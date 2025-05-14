! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> A PDE based filter

module Najafi_Yazdi_filter
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use coefs, only: coef_t
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only: pc_t, precon_factory, precon_destroy
  use bc_list, only : bc_list_t
  use neumann, only: neumann_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use pnpn_residual, only: pnpn_prs_res_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use field_registry, only: neko_field_registry
  use filter, only: filter_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_copy, field_add2
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use dofmap, only: dofmap_t
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use hsmg, only: hsmg_t
  use utils, only: neko_error, neko_type_error
  use math, only : invcol2, col2, cfill, copy
  use device_math, only: device_cfill, device_subcol3, &
                         device_cmult, device_invcol2, &
                         device_col2, device_copy
  use device, only : HOST_TO_DEVICE, device_memcpy
  use operators, only : dudxyz, div
  use projection, only : projection_t
  use time_step_controller, only : time_step_controller_t
  implicit none
  private

  ! List of all possible types created by the factory routine
  character(len=20) :: FILTER_DELTA_KNOWN_TYPES(3) = [character(len=20) :: &
       "maxGLL", &
       "averageGLL", &
       "localGLL"]

  !> A PDE based filter mapping $\rho \mapsto \tilde{\rho}$,
  !! see A. Najafi-Yazdi et al. 2015,
  !! by solving an equation
  !! of the form $\f \tilde{\rho} + alpha^2 \nabla^2 \tilde{\rho} = 
  !! \rho + beta^2 \nabla^2 \rho \f$
  type, public, extends(filter_t) :: Najafi_Yazdi_filter_t

     !> Ax
     class(ax_t), allocatable :: Ax_L, AX_R
     !> Solver results monitors ( filter )
     type(ksp_monitor_t) :: ksp_results(1)
     !> Krylov solver for the filter
     class(ksp_t), allocatable :: ksp_filt
     !> Filter Preconditioner
     class(pc_t), allocatable :: pc_filt
     !> Filter boundary conditions (they will all be Neumann, so empty)
     type(bc_list_t) :: bclst_filt
     !> b and x in equation Ax=b
     type(field_t) :: RHS, d_F_out

     ! Inputs from the user
     !> filter radius
     type(field_t) :: alpha2, beta2
     !> Transfer function at curoff wavenumber
     real(kind=rp) :: G_cutoff
     !> tolerance for PDE filter
     real(kind=rp) :: abstol_filt
     !> max iterations for PDE filter
     integer :: ksp_max_iter
     !> method for solving PDE
     character(len=:), allocatable :: ksp_solver
     !> preconditioner type
     character(len=:), allocatable :: precon_type_filt
     integer :: ksp_n, n, i
     !> If write out iteration info
     logical :: if_log
     !> Apply filter only close to elementary interfaces
     logical :: interface_only
     integer :: adjacent_idx

     !> Projection related attributes to reduce the number of iterations
     type(projection_t) :: projection
     integer :: projection_dim, projection_activ_step

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => Najafi_Yazdi_filter_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_components_field => &
          Najafi_Yazdi_filter_init_from_components_field
     procedure, pass(this) :: init_from_components_uniform => &
          Najafi_Yazdi_filter_init_from_components_uniform
     !> Destructor.
     procedure, pass(this) :: free => Najafi_Yazdi_filter_free
     !> Apply the filter
     procedure, pass(this) :: apply => Najafi_Yazdi_filter_apply
  end type Najafi_Yazdi_filter_t

contains

  !> Constructor from json.
  subroutine Najafi_Yazdi_filter_init_from_json(this, json, coef)
    class(Najafi_Yazdi_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    character(len=:), allocatable :: delta_name
    real(kind=rp) :: delta_value
    
    ! user parameters
    if (json%valid_path("filter.delta_field") .and. &
    .not. json%valid_path("filter.delta_value")) then
       call json_get(json, "filter.delta_field", delta_name)

    else if (.not. json%valid_path("filter.delta_field") .and. &
    json%valid_path("filter.delta_value")) then
       call json_get(json, "filter.delta_value", delta_value)

    else if (.not. json%valid_path("filter.delta_field") .and. &
    .not. json%valid_path("filter.delta_value")) then
       call neko_error("Please provide a delta field name or &
       &delta value to the PDE filter")
    else
       call neko_error("Please not specify a delta field name or &
       &delta value to the PDE filter together")
    end if
    
    call json_get_or_default(json, "filter.interface_only", &
         this%interface_only, .false.)
    call json_get_or_default(json, "filter.adjacent_idx", &
         this%adjacent_idx, 0)
    call json_get_or_default(json, "filter.G_cutoff", &
         this%G_cutoff, 0.5_rp)

    call json_get_or_default(json, "filter.tolerance", this%abstol_filt, &
         1.0e-10_rp)

    call json_get_or_default(json, "filter.max_iter", this%ksp_max_iter, 200)

    call json_get_or_default(json, "filter.solver", this%ksp_solver, 'cg')

    call json_get_or_default(json, "filter.preconditioner", &
         this%precon_type_filt, 'jacobi')
    
    call json_get_or_default(json, "filter.log", &
         this%if_log, .false.)
    
    call json_get_or_default(json, "filter.projection_space_size", &
         this%projection_dim, 0)
    call json_get_or_default(json, "filter.projection_hold_steps", &
         this%projection_activ_step, 0)

    call this%init_base(json, coef)
    if (json%valid_path("filter.delta_field")) then
       call Najafi_Yazdi_filter_init_from_components_field(this, &
       coef, delta_name)
    else
       call Najafi_Yazdi_filter_init_from_components_uniform(this, &
       coef, delta_value)
    end if

  end subroutine Najafi_Yazdi_filter_init_from_json

  !> Actual constructor.
  subroutine Najafi_Yazdi_filter_init_from_components_field(this, coef, &
  delta_name)
    class(Najafi_Yazdi_filter_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    character(len=:), allocatable, intent(in) :: delta_name
    real(kind=rp) :: di, dj, dk, delta_local, delta_min
    real(kind=rp) :: delta(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv)
    integer :: n, i, j, k, e, im, ip, jm, jp, km, kp
    real(kind=rp) :: delta_edge(2,2,2)
    real(kind=rp) :: pi = 4.0_rp * atan(1.0_rp)
    integer :: lx_half, ly_half, lz_half
    real(kind=rp) :: volume_element

    n = this%coef%dof%size()

    ! init the bc list (all Neuman BCs, will remain empty)
    call this%bclst_filt%init()

    ! init the b and x in Ax=b
    call this%RHS%init(this%coef%dof)
    call this%d_F_out%init(this%coef%dof)
    call this%alpha2%init(this%coef%dof)
    call this%beta2%init(this%coef%dof)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax_L, full_formulation = .false.)
    call ax_helm_factory(this%Ax_R, full_formulation = .false.)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp_filt, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol_filt)

    ! set up preconditioner
    call filter_precon_factory(this%pc_filt, this%ksp_filt, &                      
                                      this%coef, this%coef%dof, &
                                      this%coef%gs_h, &      
                                      this%bclst_filt, this%precon_type_filt)
   
    ! set up projection
    call this%projection%init(this%coef%dof%size(), this%projection_dim, &
         this%projection_activ_step)

    ! set up a delta field to generate alpha and beta
    lx_half = coef%Xh%lx / 2
    ly_half = coef%Xh%ly / 2
    lz_half = coef%Xh%lz / 2

    if (delta_name .eq. "maxGLL") then
       ! use a same length scale throughout an entire element
       ! the length scale is based on maximum GLL spacing
       do e = 1, coef%msh%nelv
          di = (coef%dof%x(lx_half, 1, 1, e) &
              - coef%dof%x(lx_half + 1, 1, 1, e))**2 &
             + (coef%dof%y(lx_half, 1, 1, e) &
              - coef%dof%y(lx_half + 1, 1, 1, e))**2 &
             + (coef%dof%z(lx_half, 1, 1, e) &
              - coef%dof%z(lx_half + 1, 1, 1, e))**2

          dj = (coef%dof%x(1, ly_half, 1, e) &
              - coef%dof%x(1, ly_half + 1, 1, e))**2 &
             + (coef%dof%y(1, ly_half, 1, e) &
              - coef%dof%y(1, ly_half + 1, 1, e))**2 &
             + (coef%dof%z(1, ly_half, 1, e) &
              - coef%dof%z(1, ly_half + 1, 1, e))**2

          dk = (coef%dof%x(1, 1, lz_half, e) &
              - coef%dof%x(1, 1, lz_half + 1, e))**2 &
             + (coef%dof%y(1, 1, lz_half, e) &
              - coef%dof%y(1, 1, lz_half + 1, e))**2 &
             + (coef%dof%z(1, 1, lz_half, e) &
              - coef%dof%z(1, 1, lz_half + 1, e))**2
          di = sqrt(di)
          dj = sqrt(dj)
          dk = sqrt(dk)
          delta(:,:,:,e) = (di * dj * dk)**(1.0_rp / 3.0_rp)
       end do
    else if (delta_name .eq. "averageGLL") then
       ! use a same length scale throughout an entire element
       ! the length scale is based on (volume)^(1/3)/(N+1)
       do e = 1, coef%msh%nelv
          volume_element = 0.0_rp
          do k = 1, coef%Xh%lx * coef%Xh%ly * coef%Xh%lz
             volume_element = volume_element + coef%B(k, 1, 1, e)
          end do
          delta(:,:,:,e) = (volume_element / (coef%Xh%lx - 1) &
            / (coef%Xh%ly - 1) / (coef%Xh%lz - 1))**(1.0_rp / 3.0_rp)
       end do
    else if (delta_name .eq. "localGLL") then
       do e = 1, coef%msh%nelv
          do k = 1, coef%Xh%lz
             km = max(1, k-1)
             kp = min(coef%Xh%lz, k+1)

             do j = 1, coef%Xh%ly
                jm = max(1, j-1)
                jp = min(coef%Xh%ly, j+1)

                do i = 1, coef%Xh%lx
                   im = max(1, i-1)
                   ip = min(coef%Xh%lx, i+1)

                   di = (coef%dof%x(ip, j, k, e) - &
                         coef%dof%x(im, j, k, e))**2 &
                      + (coef%dof%y(ip, j, k, e) - &
                         coef%dof%y(im, j, k, e))**2 &
                      + (coef%dof%z(ip, j, k, e) - &
                         coef%dof%z(im, j, k, e))**2

                   dj = (coef%dof%x(i, jp, k, e) - &
                         coef%dof%x(i, jm, k, e))**2 &
                      + (coef%dof%y(i, jp, k, e) - &
                         coef%dof%y(i, jm, k, e))**2 &
                      + (coef%dof%z(i, jp, k, e) - &
                         coef%dof%z(i, jm, k, e))**2

                   dk = (coef%dof%x(i, j, kp, e) - &
                         coef%dof%x(i, j, km, e))**2 &
                      + (coef%dof%y(i, j, kp, e) - &
                         coef%dof%y(i, j, km, e))**2 &
                      + (coef%dof%z(i, j, kp, e) - &
                         coef%dof%z(i, j, km, e))**2

                   di = sqrt(di) / (ip - im)
                   dj = sqrt(dj) / (jp - jm)
                   dk = sqrt(dk) / (kp - km)
                   delta(i,j,k,e) = (di * dj * dk)**(1.0_rp / 3.0_rp)
                end do
             end do
          end do
       end do
    else
       call neko_type_error("delta_field for filter", &
            delta_name, FILTER_DELTA_KNOWN_TYPES)
       stop
    end if

    ! set up coefficient for the laplacian on the LES and RHS
    do e = 1, this%coef%msh%nelv
       do k = 1, 2
          if (k .eq. 1) then
             km = 1
             kp = 2
          else
             km = this%coef%Xh%lz - 1
             kp = this%coef%Xh%lz
          end if
          do j = 1, 2
             if (j .eq. 1) then
                jm = 1
                jp = 2
             else
                jm = this%coef%Xh%ly - 1
                jp = this%coef%Xh%ly
             end if
             do i = 1, 2
                if (i .eq. 1) then
                   im = 1
                   ip = 2
                else
                   im = this%coef%Xh%lx - 1
                   ip = this%coef%Xh%lx
                end if

                di = (this%coef%dof%x(ip, j, k, e) - &
                   this%coef%dof%x(im, j, k, e))**2 &
                   + (this%coef%dof%y(ip, j, k, e) - &
                   this%coef%dof%y(im, j, k, e))**2 &
                   + (this%coef%dof%z(ip, j, k, e) - &
                   this%coef%dof%z(im, j, k, e))**2

                dj = (this%coef%dof%x(i, jp, k, e) - &
                   this%coef%dof%x(i, jm, k, e))**2 &
                   + (this%coef%dof%y(i, jp, k, e) - &
                   this%coef%dof%y(i, jm, k, e))**2 &
                   + (this%coef%dof%z(i, jp, k, e) - &
                   this%coef%dof%z(i, jm, k, e))**2

                dk = (this%coef%dof%x(i, j, kp, e) - &
                   this%coef%dof%x(i, j, km, e))**2 &
                   + (this%coef%dof%y(i, j, kp, e) - &
                   this%coef%dof%y(i, j, km, e))**2 &
                   + (this%coef%dof%z(i, j, kp, e) - &
                   this%coef%dof%z(i, j, km, e))**2
                di = sqrt(di)
                dj = sqrt(dj)
                dk = sqrt(dk)
                delta_edge(i,j,k) = (di * dj * dk)**(1.0_rp / 3.0_rp)
             end do
          end do
       end do
       delta_min = minval(delta_edge)
       this%beta2%x(:,:,:,e) = - delta_min * delta_min /pi/pi

       do k = 1, this%coef%Xh%lz
          do j = 1, this%coef%Xh%ly
             do i = 1, this%coef%Xh%lz
                delta_local = delta(i,j,k,e)
                this%alpha2%x(i,j,k,e) = -1.0_rp * delta_local * delta_local &
                     / pi / pi * (1.0_rp - 1.0_rp/this%G_cutoff * &
                     (1.0_rp - delta_min * delta_min / &
                     delta_local / delta_local))
             end do
          end do
       end do
       
    end do
    if (this%interface_only) then
       this%beta2%x(2 + this%adjacent_idx:this%coef%Xh%lx-1-this%adjacent_idx, &
         2 + this%adjacent_idx : this%coef%Xh%ly - 1 - this%adjacent_idx, &
         2 + this%adjacent_idx : this%coef%Xh%lz - 1 - this%adjacent_idx, &
         :) = 0.0_rp
       this%alpha2%x(2 + this%adjacent_idx:this%coef%Xh%lx-1-this%adjacent_idx, &
         2 + this%adjacent_idx : this%coef%Xh%ly - 1 - this%adjacent_idx, &
         2 + this%adjacent_idx : this%coef%Xh%lz - 1 - this%adjacent_idx, &
         :) = 0.0_rp
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%beta2%x, this%beta2%x_d, this%beta2%dof%size(), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%alpha2%x, this%alpha2%x_d, this%alpha2%dof%size(), &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine Najafi_Yazdi_filter_init_from_components_field

  subroutine Najafi_Yazdi_filter_init_from_components_uniform(this, coef, &
  delta_value)
    class(Najafi_Yazdi_filter_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: delta_value
    real(kind=rp) :: di, dj, dk, delta_min
    integer :: n, i, j, k, e, im, ip, jm, jp, km, kp
    real(kind=rp) :: delta_edge(2,2,2)
    real(kind=rp) :: pi = 4.0_rp * atan(1.0_rp)

    n = this%coef%dof%size()

    ! init the bc list (all Neuman BCs, will remain empty)
    call this%bclst_filt%init()

    ! init the b and x in Ax=b
    call this%RHS%init(this%coef%dof)
    call this%d_F_out%init(this%coef%dof)
    call this%alpha2%init(this%coef%dof)
    call this%beta2%init(this%coef%dof)

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax_L, full_formulation = .false.)
    call ax_helm_factory(this%Ax_R, full_formulation = .false.)

    ! set up krylov solver
    call krylov_solver_factory(this%ksp_filt, n, this%ksp_solver, &
         this%ksp_max_iter, this%abstol_filt)

    ! set up preconditioner
    call filter_precon_factory(this%pc_filt, this%ksp_filt, &                      
                                      this%coef, this%coef%dof, &
                                      this%coef%gs_h, &      
                                      this%bclst_filt, this%precon_type_filt)
   
    ! set up projection
    call this%projection%init(this%coef%dof%size(), this%projection_dim, &
         this%projection_activ_step)

    ! set up coefficient for the laplacian on the LES and RHS
    do e = 1, this%coef%msh%nelv
       do k = 1, 2
          if (k .eq. 1) then
             km = 1
             kp = 2
          else
             km = this%coef%Xh%lz - 1
             kp = this%coef%Xh%lz
          end if
          do j = 1, 2
             if (j .eq. 1) then
                jm = 1
                jp = 2
             else
                jm = this%coef%Xh%ly - 1
                jp = this%coef%Xh%ly
             end if
             do i = 1, 2
                if (i .eq. 1) then
                   im = 1
                   ip = 2
                else
                   im = this%coef%Xh%lx - 1
                   ip = this%coef%Xh%lx
                end if

                di = (this%coef%dof%x(ip, j, k, e) - &
                   this%coef%dof%x(im, j, k, e))**2 &
                   + (this%coef%dof%y(ip, j, k, e) - &
                   this%coef%dof%y(im, j, k, e))**2 &
                   + (this%coef%dof%z(ip, j, k, e) - &
                   this%coef%dof%z(im, j, k, e))**2

                dj = (this%coef%dof%x(i, jp, k, e) - &
                   this%coef%dof%x(i, jm, k, e))**2 &
                   + (this%coef%dof%y(i, jp, k, e) - &
                   this%coef%dof%y(i, jm, k, e))**2 &
                   + (this%coef%dof%z(i, jp, k, e) - &
                   this%coef%dof%z(i, jm, k, e))**2

                dk = (this%coef%dof%x(i, j, kp, e) - &
                   this%coef%dof%x(i, j, km, e))**2 &
                   + (this%coef%dof%y(i, j, kp, e) - &
                   this%coef%dof%y(i, j, km, e))**2 &
                   + (this%coef%dof%z(i, j, kp, e) - &
                   this%coef%dof%z(i, j, km, e))**2
                di = sqrt(di)
                dj = sqrt(dj)
                dk = sqrt(dk)
                delta_edge(i,j,k) = (di * dj * dk)**(1.0_rp / 3.0_rp)
             end do
          end do
       end do
       delta_min = minval(delta_edge)
       this%beta2%x(:,:,:,e) = - delta_min * delta_min /pi/pi

       do k = 1, this%coef%Xh%lz
          do j = 1, this%coef%Xh%ly
             do i = 1, this%coef%Xh%lz
                this%alpha2%x(i,j,k,e) = -1.0_rp * delta_value * delta_value &
                     / pi / pi * (1.0_rp - 1.0_rp/this%G_cutoff * &
                     (1.0_rp - delta_min * delta_min / &
                     delta_value / delta_value))
             end do
          end do
       end do
       
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%beta2%x, this%beta2%x_d, this%beta2%dof%size(), &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%alpha2%x, this%alpha2%x_d, this%alpha2%dof%size(), &
            HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine Najafi_Yazdi_filter_init_from_components_uniform

  !> Destructor.
  subroutine Najafi_Yazdi_filter_free(this)
    class(Najafi_Yazdi_filter_t), intent(inout) :: this

    if (allocated(this%Ax_L)) then
       deallocate(this%Ax_L)
    end if

    if (allocated(this%Ax_R)) then
       deallocate(this%Ax_R)
    end if                                                                                         
                                                                                
    if (allocated(this%pc_filt)) then                                                
       call precon_destroy(this%pc_filt)                                             
       deallocate(this%pc_filt)                                                      
    end if                    

    call this%bclst_filt%free()
    call this%projection%free()
    call this%RHS%free()
    call this%d_F_out%free()
    call this%alpha2%free()
    call this%beta2%free()
    call this%ksp_filt%free()

    call this%free_base()

  end subroutine Najafi_Yazdi_filter_free

  !> Apply the filter
  !! @param F_out filtered field
  !! @param F_in unfiltered field
  subroutine Najafi_Yazdi_filter_apply(this, F_out, F_in, tstep, dt_controller)
    class(Najafi_Yazdi_filter_t), intent(inout) :: this
    type(field_t), intent(in) :: F_in
    type(field_t), intent(inout) :: F_out
    integer, intent(in), optional :: tstep
    type(time_step_controller_t), intent(in), optional :: dt_controller
    integer :: n, i
    character(len=LOG_SIZE) :: log_buf

    n = this%coef%dof%size()

    ! set up Helmholtz operators for \rho + beta^2 \nabla^2 \rho
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%coef%h1_d, this%beta2%x_d, n)
       call device_cfill(this%coef%h2_d, 1.0_rp, n)
    else
       ! h1 is already negative in its definition
       call copy(this%coef%h1, this%beta2%x, n)
       ! ax_helm includes the mass matrix in h2
       call cfill(this%coef%h2, 1.0_rp, n)
    end if
    this%coef%ifh2 = .true.

    ! compute the A(F_in) component of the RHS 
    ! use F_out as a temporal array
    ! First, construct \rho + beta^2 \nabla^2 \rho
    call this%Ax_R%compute(F_out%x, F_in%x, this%coef, this%coef%msh, &
        this%coef%Xh)
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol2(F_out%x_d, this%coef%B_d, n)
    else
       call invcol2(F_out%x, this%coef%B, n)
    end if

    ! gather scatter
    call this%coef%gs_h%op(F_out, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(F_out%x_d, this%coef%mult_d, n)
    else
       call col2(F_out%x, this%coef%mult, n)
    end if

    ! set up Helmholtz operators for equation solving
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%coef%h1_d, this%alpha2%x_d, n)
      !  call device_cfill(this%coef%h2_d, 1.0_rp, n)
    else
       ! h1 is already negative in its definition
       call copy(this%coef%h1, this%alpha2%x, n)
      !  ! ax_helm includes the mass matrix in h2
      !  call cfill(this%coef%h2, 1.0_rp, n)
    end if
    this%coef%ifh2 = .true.

    call this%Ax_L%compute(this%RHS%x, F_out%x, this%coef, this%coef%msh, &
        this%coef%Xh)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_subcol3(this%RHS%x_d, F_out%x_d, this%coef%B_d, n)
       call device_cmult(this%RHS%x_d, -1.0_rp, n)
    else
       do i = 1, n
          ! mass matrix should be included here
          this%RHS%x(i,1,1,1) = F_out%x(i,1,1,1) * this%coef%B(i,1,1,1) &
              - this%RHS%x(i,1,1,1)
       end do
    end if

    ! gather scatter
    call this%coef%gs_h%op(this%RHS, GS_OP_ADD)

    ! set BCs
    call this%bclst_filt%apply_scalar(this%RHS%x, n)

    call this%projection%pre_solving(this%RHS%x, tstep, this%coef, n, &
         dt_controller)

    ! Solve Helmholtz equation
    call profiler_start_region('filter solve')
    this%ksp_results(1) = &
         this%ksp_filt%solve(this%Ax_L, this%d_F_out, this%RHS%x, n, &
         this%coef, this%bclst_filt, this%coef%gs_h)

    call profiler_end_region

    call this%projection%post_solving(this%d_F_out%x, this%Ax_L, this%coef, &
         this%bclst_filt, this%coef%gs_h, n, tstep, dt_controller)

    ! add result
    call field_add2(F_out, this%d_F_out)
    ! update preconditioner (needed?)
    call this%pc_filt%update()
    
    if (this%if_log) then
       ! write it all out
       call neko_log%message('Filter')

       write(log_buf, '(A,A,A)') 'Iterations:   ',&
             'Start residual:     ', 'Final residual:'
       call neko_log%message(log_buf)
       write(log_buf, '(I11,3x, E15.7,5x, E15.7)') this%ksp_results%iter, &
             this%ksp_results%res_start, this%ksp_results%res_final
       call neko_log%message(log_buf)
    end if

  end subroutine Najafi_Yazdi_filter_apply

  !> Initialize a Krylov preconditioner
  subroutine filter_precon_factory(pc, ksp, coef, dof, gs, bclst, &
       pctype)
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type (pcp => pc)
      type is (jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
      type is (hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst, &
                  trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    end select

    call ksp%set_pc(pc)

  end subroutine filter_precon_factory

end module Najafi_Yazdi_filter