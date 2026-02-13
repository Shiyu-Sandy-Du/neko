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
!> A simulation component that computes entropy_viscosity
!! The values are stored in the field registry under the name 'entropy_viscosity'

module entropy_viscosity
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use json_utils, only : json_get, json_get_or_default
  use field, only : field_t, field_ptr_t
  use field_series, only : field_series_t
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t
  use advection, only : advection_t
  use time_scheme_controller, only : time_scheme_controller_t
  use coefs, only : coef_t
  use field_series, only : field_series_ptr_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use fluid_pnpn, only : fluid_pnpn_t
  use scalars, only : scalars_t
  use utils, only : neko_error
  use elementwise_filter, only : elementwise_filter_t
  use field_math, only : field_col3, field_copy, field_absval, field_rzero, &
                         field_cmult, field_sub2, field_col2, field_cadd2, &
                         field_invcol2, field_sqrt, field_add2, field_addcol3, &
                         field_invcol2_nonzero, field_add3, field_pwmin2
  use math, only : NEKO_EPS, invcol2, col2, glsum, glsc2, glmax, glmin
  use device_math, only : device_invcol2, device_col2, device_glsum, &
                          device_glsc2, device_glmax, device_glmin
  use tensor, only : dottnsr_3d, dot1tnsr_3d, maxnorm_3d
  use gather_scatter, only : GS_OP_ADD
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: entropy_viscosity_t
     !> coefficient
     real(kind=rp) :: c_E
     !> Upper bound coefficient
     real(kind=rp) :: c_max = 0.5_rp
     type(field_t) :: h_k 
     !> A low pass filter for the field
     type(elementwise_filter_t) :: filter
     logical :: if_filter = .false.
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> work array for the temporal derivative
     type(field_t), allocatable :: wa(:)
     !> Scalar field
     integer :: n_scalars = 0
     type(field_ptr_t), allocatable :: s(:)
     type(field_t), allocatable :: E(:)
     type(field_series_t), allocatable :: Elag(:)
     !> coef
     type(coef_t), pointer :: coef
     !> Some field to be used
     type(field_t) :: h2
     real(kind=rp) :: volume_domain

     !> Scaling option
     character(len=:), allocatable :: scaling_option

     !> Tolerance coefficient for a minimum entropy fluctuation level
     !! with respect to the maximum fluctuation in the whole domain
     real(kind=rp) :: tol_coef

     !> entropy_viscosity fields.
     type(field_ptr_t), allocatable :: entropy_viscosity(:)

     !> Residual.
     type(field_ptr_t), allocatable :: D(:)

     !> Output writer.
     type(field_writer_t) :: writer

     !> A pointer pointing to the time scheme
     type(fluid_pnpn_t), pointer :: fluid
     type(scalars_t), pointer :: scalars
     class(rhs_maker_bdf_t), pointer :: makebdf
     class(rhs_maker_ext_t), pointer :: makeext
     class(advection_t), pointer :: adv
     type(time_scheme_controller_t), pointer :: ext_bdf

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => entropy_viscosity_init_from_json
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => entropy_viscosity_init_common
     !> Destructor.
     procedure, pass(this) :: free => entropy_viscosity_free
     !> Part of the residual viscosity computation before the time stepping
     procedure, pass(this) :: preprocess_ => entropy_viscosity_preprocess
     !> Part of the residual viscosity computation after the time stepping
     procedure, pass(this) :: compute_ => entropy_viscosity_compute
  end type entropy_viscosity_t

contains

  !> Constructor from json.
  subroutine entropy_viscosity_init_from_json(this, json, case)
    class(entropy_viscosity_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=:), allocatable :: scaling_option

    call this%init_base(json, case)
    
    call json_get(json, "c_E", this%c_E)
    call json_get(json, "scaling_option", scaling_option)
    this%scaling_option = trim(scaling_option)

    select case (trim(scaling_option))
    case ("global_average")
    case ("global_minmax")
   !  case ("elementwise")
   !     call json_get_or_default(json, "tol_coef", this%tol_coef, 1e-3_rp)
    case default
       call neko_error("Invalide scaling option for entropy viscosity")
    end select

    

    call this%init_common(json, case)
  end subroutine entropy_viscosity_init_from_json

  !> Common part of constructors.
  subroutine entropy_viscosity_init_common(this, json, case)
    class(entropy_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20), allocatable :: fields(:)
    integer :: e, k
    real(kind=rp) :: volume_element, diameter

    if (allocated(case%scalars)) then
       this%scalars => case%scalars
       this%n_scalars = size(this%scalars%scalar_fields)
    else
       this%n_scalars = 0
    end if
    allocate(fields(2*(1+this%n_scalars)))
    fields(1) = 'entr_visc_vel'
    fields(2) = 'entr_res_vel'
    do k = 1, this%n_scalars
       write(fields(2*k+1), '(A,I0)') 'entr_visc_s', k
       write(fields(2*k+2), '(A,I0)') 'entr_res_s', k
    end do
    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/entropy_viscosity.f90\the registry.
    call json%add("fields", fields)
    call this%writer%init(json, case)

    this%coef => case%fluid%c_Xh
    
    ! Set up the filter
    if (json%valid_path("filter")) then
       this%if_filter = .true.
       call this%filter%init(json, this%coef)
       this%filter%transfer(this%coef%dof%xh%lx) = 0.0_rp ! filter out the highest order mode
       call this%filter%build_1d()
    end if

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
      this%makebdf => f1%makebdf
      this%makeext => f1%makeabf
      this%adv => f1%adv
      this%ext_bdf => f1%ext_bdf
    class default
      call neko_error("For fluid, entropy &
      &viscosity currently only support pnpn scheme")
    end select

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")


    call this%h2%init(this%u%dof)
    call this%h_k%init(this%u%dof)

    if (this%n_scalars .ne. 0) then
       allocate(this%s(this%n_scalars))
    end if
    
    allocate(this%E(1+this%n_scalars))
    allocate(this%Elag(1+this%n_scalars))
    allocate(this%wa(1+this%n_scalars))

    allocate(this%D(1+this%n_scalars))
    
    allocate(this%entropy_viscosity(1+this%n_scalars))

    do k = 1, 1+this%n_scalars
       this%entropy_viscosity(k)%ptr => &
              neko_field_registry%get_field(fields(2*k-1))
       this%D(k)%ptr => &
              neko_field_registry%get_field(fields(2*k))

       call this%E(k)%init(this%u%dof)
       call this%Elag(k)%init(this%E(k), 2)
       call this%wa(k)%init(this%u%dof)
       
       if (k .le. this%n_scalars) then
          this%s(k)%ptr => this%scalars%scalar_fields(k)%s
       end if
    end do

    do e = 1, this%coef%msh%nelv
      !  volume_element = 0.0_rp
      !  do k = 1, this%coef%Xh%lx * this%coef%Xh%ly * this%coef%Xh%lz
      !     volume_element = volume_element + this%coef%B(k, 1, 1, e)
      !  end do
      !  this%h2%x(:,:,:,e) = volume_element**(1.0_rp/3.0_rp) * &
      !                       volume_element**(1.0_rp/3.0_rp) / &
      !                       (this%coef%Xh%lx-1.0_rp) / &
      !                       (this%coef%Xh%lx-1.0_rp)
       diameter = this%coef%msh%elements(e)%e%diameter()
       this%h2%x(:,:,:,e) = diameter * diameter / &
                            (this%coef%Xh%lx-1.0_rp) / &
                            (this%coef%Xh%lx-1.0_rp)
       this%h_k%x(:,:,:,e) = this%coef%msh%elements(e)%e%diameter()
    end do
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%h2%x, this%h2%x_d, this%u%dof%size(), &
                             HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%h_k%x, this%h_k%x_d, this%u%dof%size(), &
                             HOST_TO_DEVICE, sync = .false.)
    end if
   
    if (this%scaling_option .eq. "global_average") then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          this%volume_domain = device_glsum(this%coef%B_d, this%u%dof%size())
       else
          this%volume_domain = glsum(this%coef%B, this%u%dof%size())
       end if
    end if

  end subroutine entropy_viscosity_init_common

  !> Destructor.
  subroutine entropy_viscosity_free(this)
    class(entropy_viscosity_t), intent(inout) :: this
    call this%free_base()
  end subroutine entropy_viscosity_free

  !> Part of the entropy_viscosity computation before the time stepping.
  !! @param time The time state.
  subroutine entropy_viscosity_preprocess(this, time)
    class(entropy_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, n

    ! Time lag part for the BDF scheme of dE/dt
    associate(wa_vel => this%wa(1), E_vel => this%E(1), &
              Elag_vel => this%Elag(1), &
              coef => this%coef, &
              rho => this%fluid%rho, dt => time%dt, &
              makebdf => this%makebdf, &
              ext_bdf => this%ext_bdf)

     n = wa_vel%dof%size()
     call field_rzero(wa_vel)
     call makebdf%compute_scalar(Elag_vel, wa_vel%x, E_vel, coef%B, &
               rho%x(1,1,1,1), dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
     if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_invcol2(wa_vel%x_d, coef%B_d, n)
     else
        call invcol2(wa_vel%x, coef%B, n)
     end if
     call Elag_vel%update()
     
    end associate
    
    do i = 1, this%n_scalars
       ! Time lag part for the BDF scheme of dE/dt
       associate(wa_s_i => this%wa(i+1), &
                 coef => this%coef, &
                 rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
                 makebdf => this%makebdf, E_s_i => this%E(i+1), &
                 Elag_s_i => this%Elag(i+1), ext_bdf => this%ext_bdf)
      
       call field_rzero(wa_s_i)
       call makebdf%compute_scalar(Elag_s_i, wa_s_i%x, E_s_i, coef%B, rho%x(1,1,1,1), &
               dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(wa_s_i%x_d, coef%B_d, n)
       else
          call invcol2(wa_s_i%x, coef%B, n)
       end if
       call Elag_s_i%update()

       end associate
    end do

  end subroutine entropy_viscosity_preprocess

  !> Part of the entropy_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine entropy_viscosity_compute(this, time)
    class(entropy_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ta ! temporal array
    type(field_ptr_t) :: fs(this%n_scalars)
    type(field_t), pointer :: fu, fv, fw

    ! Upper bound
    type(field_t), pointer :: local_max_ws ! ws for wave speed

    ! Variables to be used if the scaling option is global average
    type(field_t), pointer :: E_vel_var
    type(field_ptr_t) :: E_s_var(this%n_scalars)
    real(kind=rp) :: E_vel_avg
    real(kind=rp) :: E_s_avg
    integer :: E_vel_var_index
    integer :: E_s_var_indices(this%n_scalars)

    ! Variables to be used if the scaling option is elementwise
    type(field_t), pointer :: E_vel_avg_field
    type(field_ptr_t) :: E_s_avg_field(this%n_scalars)
    type(field_t), pointer :: B_elem
    integer :: E_vel_avg_index
    integer :: E_s_avg_indices(this%n_scalars)
    integer :: B_elem_index
    real(kind=rp) :: tol

    integer :: temp_indices(2)
    integer :: filt_field_indices(3+this%n_scalars)

    integer :: i, j, n
    real(kind=rp) :: scaling_factor

    call neko_scratch_registry%request_field(ta, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(local_max_ws, &
         temp_indices(2), .false.)

    if (this%scaling_option .eq. "global_average") then
       call neko_scratch_registry%request_field(E_vel_var, &
                                             E_vel_var_index, .false.)
       do i = 1, this%n_scalars
          call neko_scratch_registry%request_field(E_s_var(i)%ptr, &
                                             E_s_var_indices(i), .false.)
       end do

   !  else if (this%scaling_option .eq. "elementwise") then
   !     call neko_scratch_registry%request_field(E_vel_avg_field, &
   !                                              E_vel_avg_index, .false.)

   !     do i = 1, this%n_scalars
   !        call neko_scratch_registry%request_field(E_s_avg_field(i)%ptr, &
   !                                                 E_s_avg_indices(i), .false.)
   !     end do
   !     call neko_scratch_registry%request_field(B_elem, B_elem_index, .false.)
    end if

    call neko_scratch_registry%request_field(fu, filt_field_indices(1), .false.)
    call neko_scratch_registry%request_field(fv, filt_field_indices(2), .false.)
    call neko_scratch_registry%request_field(fw, filt_field_indices(3), .false.)
    do i = 1, this%n_scalars
       call neko_scratch_registry%request_field(fs(i)%ptr, &
            filt_field_indices(3+i), .false.)
    end do

    ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
    associate(u => this%u, v => this%v, w => this%w, E_vel => this%E(1), &
             ext_bdf => this%ext_bdf, &
             dt => time%dt, coef => this%coef, wa_vel => this%wa(1), &
             D_vel => this%D(1)%ptr, gs => this%coef%gs_h, &
             adv => this%adv, &
             Xh => this%coef%Xh, &
             entropy_viscosity_vel => this%entropy_viscosity(1)%ptr)

    n = u%dof%size()

    if (this%if_filter) then
      ! filter the velocity magnitude all together
      call field_rzero(ta)
      call field_addcol3(ta, u, u)
      call field_addcol3(ta, v, v)
      call field_addcol3(ta, w, w)
      call field_sqrt(ta)
      call maxnorm_3d(local_max_ws%x, ta%x, coef%Xh%lx, coef%msh%nelv)
      call this%filter%apply(fu, ta)
      call field_sub2(fu, ta)
      call gs%op(fu, GS_OP_ADD)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(fu%x_d, coef%mult_d, n)
      else
         call col2(fu%x, coef%mult, n)
      end if
      call field_copy(E_vel, fu)
      call field_col2(E_vel, E_vel)
    else
      call field_col3(ta, u, u)
      call field_copy(E_vel, ta)
      call field_col3(ta, v, v)
      call field_add2(E_vel, ta)
      call field_col3(ta, w, w)
      call field_add2(E_vel, ta)
      call field_copy(ta, E_vel)
      call field_sqrt(ta)
      call maxnorm_3d(local_max_ws%x, ta%x, coef%Xh%lx, coef%msh%nelv)
    end if

    call field_col2(local_max_ws, this%h_k)
    ! now local_max_ws is a work array for the upper bound of the viscosity
    call field_cmult(local_max_ws, this%c_max)

    call field_copy(ta, E_vel)
    call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
    call field_sub2(ta, wa_vel)
    call field_copy(D_vel, ta)

    ! advection part
    call field_rzero(ta)
    call adv%compute_scalar(u, v, w, E_vel, ta, &
         Xh, coef, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_invcol2(ta%x_d, coef%B_d, n)
    else
       call invcol2(ta%x, coef%B, n)
    end if
    call gs%op(ta, GS_OP_ADD)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(ta%x_d, coef%mult_d, n)
    else
       call col2(ta%x, coef%mult, n)
    end if
    call field_sub2(D_vel, ta, n)
    call field_copy(entropy_viscosity_vel, D_vel)
    call field_absval(entropy_viscosity_vel)

    if (this%scaling_option .eq. "global_average") then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          E_vel_avg = - device_glsc2(E_vel%x_d, coef%B_d, n) / this%volume_domain
       else
          E_vel_avg = - glsc2(E_vel%x, coef%B, n) / this%volume_domain
       end if
       call field_cadd2(E_vel_var, E_vel, E_vel_avg)
       call field_absval(E_vel_var)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          scaling_factor = this%c_E / device_glmax(E_vel_var%x_d, n)
       else
          scaling_factor = this%c_E / glmax(E_vel_var%x, n)
       end if

       call field_cmult(entropy_viscosity_vel, &
         scaling_factor)

    else if (this%scaling_option .eq. "global_minmax") then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          scaling_factor = device_glmax(E_vel%x_d, n) - &
                           device_glmin(E_vel%x_d, n)
       else
          scaling_factor = glmax(E_vel%x, n) - glmin(E_vel%x, n)
       end if
       if (scaling_factor .lt. NEKO_EPS) then
          scaling_factor = 0.0_rp
       else
          scaling_factor = this%c_E / scaling_factor
       end if
       call field_cmult(entropy_viscosity_vel, scaling_factor)
   !  ! Temporal implementation of elementwise scaling
   !  else if (this%scaling_option .eq. "elementwise") then
   !     call dottnsr_3d(E_vel_avg_field%x, E_vel%x, &
   !             coef%B, coef%Xh%lx, coef%msh%nelv)
   !     call dot1tnsr_3d(B_elem%x, coef%B, coef%Xh%lx, coef%msh%nelv)
   !     call field_invcol2(E_vel_avg_field, B_elem)
   !     call field_cmult(E_vel_avg_field, -1.0_rp)

   !     call field_add3(E_vel_var, E_vel, E_vel_avg_field)
   !     call field_absval(E_vel_var)

   !     call field_cmult(entropy_viscosity_vel, &
   !             this%c_E)

   !     call maxnorm_3d(ta%x, E_vel_var%x, coef%Xh%lx, coef%msh%nelv)
   !     if (NEKO_BCKND_DEVICE .eq. 1) then
   !        tol = this%tol_coef * (device_glmax(E_vel%x_d, E_vel%dof%size()) - &
   !                               device_glmin(E_vel%x_d, E_vel%dof%size()))
   !     else
   !        tol = this%tol_coef * (glmax(E_vel%x, E_vel%dof%size()) - &
   !                               glmin(E_vel%x, E_vel%dof%size()))
   !     end if 
   !     write(*,*) "xxx", tol
   !     call field_invcol2_nonzero(entropy_viscosity_vel, ta, tol)
    end if


    call field_col2(entropy_viscosity_vel, this%h2)
    call field_pwmin2(entropy_viscosity_vel, local_max_ws)

   end associate

    do i = 1, this%n_scalars
       ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
       associate(s_i => this%s(i)%ptr, &
                 fs_i => fs(i)%ptr, E_s_i => this%E(i+1), &
                 ext_bdf => this%ext_bdf, &
                 dt => time%dt, coef => this%coef, wa_s_i => this%wa(i+1), &
                 D_s_i => this%D(i+1)%ptr, gs => this%coef%gs_h, &
                 adv => this%adv, &
                 u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
                 E_s_var_i => E_s_var(i)%ptr, &
                 E_s_avg_field_i => E_s_avg_field(i)%ptr, &
                 entropy_viscosity_i => this%entropy_viscosity(i+1)%ptr)

       n = s_i%dof%size()

       if (this%if_filter) then
         call this%filter%apply(fs_i, s_i)
         call field_sub2(fs_i, s_i)
         call gs%op(fs_i, GS_OP_ADD)
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_col2(fs_i%x_d, coef%mult_d, n)
         else
            call col2(fs_i%x, coef%mult, n)
         end if
         call field_copy(E_s_i, fs_i)
       else
         call field_copy(E_s_i, s_i)
       end if

       ! Take the square as the entropy
       call field_col2(E_s_i, E_s_i)
       call field_copy(ta, E_s_i)
       call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
       call field_sub2(ta, wa_s_i)
       call field_copy(D_s_i, ta)

       ! advection part
       call field_rzero(ta, n)
       call adv%compute_scalar(u, v, w, E_s_i, ta, &
              Xh, coef, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(ta%x_d, coef%B_d, n)
       else
          call invcol2(ta%x, coef%B, n)
       end if
       call gs%op(ta, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(ta%x_d, coef%mult_d, n)
       else
          call col2(ta%x, coef%mult, n)
       end if
       call field_sub2(D_s_i, ta, n)
       call field_copy(entropy_viscosity_i, D_s_i)
       call field_absval(entropy_viscosity_i)
       
       if (this%scaling_option .eq. "global_average") then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             E_s_avg = - device_glsc2(E_s_i%x_d, coef%B_d, n) / this%volume_domain
          else
             E_s_avg = - glsc2(E_s_i%x, coef%B, n) / this%volume_domain
          end if
          call field_cadd2(E_s_var_i, E_s_i, E_s_avg)
          call field_absval(E_s_var_i)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             scaling_factor = this%c_E / device_glmax(E_s_var_i%x_d, n)
          else
             scaling_factor = this%c_E / glmax(E_s_var_i%x, n)
          end if

          call field_cmult(entropy_viscosity_i, &
               scaling_factor)

       else if (this%scaling_option .eq. "global_minmax") then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             scaling_factor = device_glmax(E_s_i%x_d, n) - &
                              device_glmin(E_s_i%x_d, n)
          else
             scaling_factor = glmax(E_s_i%x, n) - glmin(E_s_i%x, n)
          end if

          if (scaling_factor .lt. NEKO_EPS) then
             scaling_factor = 0.0_rp
          else
             scaling_factor = this%c_E / scaling_factor
          end if

          call field_cmult(entropy_viscosity_i, scaling_factor)
      !  else if (this%scaling_option .eq. "elementwise") then
         !  call dottnsr_3d(E_s_avg_field_i%x, E_s_i%x, &
         !       coef%B, coef%Xh%lx, coef%msh%nelv)
         !  call dot1tnsr_3d(B_elem%x, coef%B, coef%Xh%lx, coef%msh%nelv)
         !  call field_invcol2(E_s_avg_field_i, B_elem)
         !  call field_cmult(E_s_avg_field_i, -1.0_rp)

         !  call field_add3(E_s_var_i, E_s_i, E_s_avg_field_i)
         !  call field_absval(E_s_var_i)

         !  call field_cmult(entropy_viscosity_i, &
         !       this%c_E)

         !  call maxnorm_3d(ta%x, E_s_var_i%x, coef%Xh%lx, coef%msh%nelv)
         !  if (NEKO_BCKND_DEVICE .eq. 1) then
         !     tol = this%tol_coef * (device_glmax(E_s_i%x_d, E_s_i%dof%size()) - &
         !                        device_glmin(E_s_i%x_d, E_s_i%dof%size()))
         !  else
         !     tol = this%tol_coef * (glmax(E_s_i%x, E_s_i%dof%size()) - &
         !                        glmin(E_s_i%x, E_s_i%dof%size()))
         !  end if 
         !  call field_invcol2_nonzero(entropy_viscosity_i, ta, tol)

       end if

       call field_col2(entropy_viscosity_i, this%h2)
       call field_pwmin2(entropy_viscosity_i, local_max_ws)

       end associate
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
    call neko_scratch_registry%relinquish_field(filt_field_indices)

    if (this%scaling_option .eq. "global_average") then
       call neko_scratch_registry%relinquish_field(E_vel_var_index)
       call neko_scratch_registry%relinquish_field(E_s_var_indices)
   !  else if (this%scaling_option .eq. "elementwise") then
   !     call neko_scratch_registry%relinquish_field(E_vel_avg_index)
   !     call neko_scratch_registry%relinquish_field(E_s_avg_indices)
   !     call neko_scratch_registry%relinquish_field(B_elem_index)
    end if

  end subroutine entropy_viscosity_compute

end module entropy_viscosity
