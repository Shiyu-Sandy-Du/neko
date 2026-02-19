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
     real(kind=rp) :: c_max
     type(field_t) :: h_k 
     !> The power coefficient for the elementwise filter
     real(kind=rp) :: power_coef = 1.5_rp
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

     !> Tolerance coefficient for a minimum entropy fluctuation level
     !! with respect to the maximum fluctuation in the whole domain
     real(kind=rp) :: tol_coef

     !> entropy_viscosity fields.
     type(field_ptr_t), allocatable :: entropy_viscosity(:)

     !> Residual.
     type(field_ptr_t), allocatable :: D(:)

     !> Upper bound
     type(field_ptr_t), allocatable :: ev_cap ! ws for wave speed

     !> Output writer.
     type(field_writer_t) :: writer

     !> A pointer pointing to the time scheme
     type(fluid_pnpn_t), pointer :: fluid
     type(scalars_t), pointer :: scalars
     class(rhs_maker_bdf_t), pointer :: makebdf
     class(rhs_maker_ext_t), pointer :: makeext
     ! adv should be forced to be no-dealised version for nodal operation
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

    call this%init_base(json, case)
    
    call json_get(json, "c_E", this%c_E)
    call json_get_or_default(json, "tol_coef", this%tol_coef, 1e-2_rp)
    call json_get_or_default(json, "c_max", this%c_max, 0.5_rp)

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
    allocate(fields(1+2*(1+this%n_scalars)))
    fields(1) = 'entr_visc_vel'
    fields(2) = 'entr_res_vel'
    do k = 1, this%n_scalars
       write(fields(2*k+1), '(A,I0)') 'entr_visc_s', k
       write(fields(2*k+2), '(A,I0)') 'entr_res_s', k
    end do
    fields(1+2*(1+this%n_scalars)) = 'max_entr_visc'

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/entropy_viscosity.f90\the registry.
    call json%add("fields", fields)
    call this%writer%init(json, case)

    this%coef => case%fluid%c_Xh
    
    ! Set up the filter
    if (json%valid_path("filter")) then
       this%if_filter = .true.
       call this%filter%init(json, this%coef)
      !  this%filter%transfer(this%coef%dof%xh%lx) = 0.0_rp ! filter out the highest order mode
       ! give the weight of around 0.2 to the second highest mode while keeping the kernel smooth
       do k = 1, this%coef%Xh%lx
          this%filter%transfer(k) = ((k - 1.0_rp) / (this%coef%Xh%lx - 1.0_rp)) &
                                 ** ((this%coef%Xh%lx - 1.0_rp) * this%power_coef)
          this%filter%transfer(k) = 1.0_rp - this%filter%transfer(k)
       end do
       call this%filter%build_1d()
    end if

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
      this%makebdf => f1%makebdf
      this%makeext => f1%makeabf
      this%ext_bdf => f1%ext_bdf
      this%adv => f1%adv
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
    allocate(this%ev_cap)

    do k = 1, 1+this%n_scalars
       this%entropy_viscosity(k)%ptr => &
              neko_field_registry%get_field(fields(2*k-1))
       this%D(k)%ptr => &
              neko_field_registry%get_field(fields(2*k))

       call field_rzero(this%entropy_viscosity(k)%ptr)
       call field_rzero(this%D(k)%ptr)

       call this%E(k)%init(this%u%dof)
       call this%Elag(k)%init(this%E(k), 2)
       call this%wa(k)%init(this%u%dof)
       
       if (k .le. this%n_scalars) then
          this%s(k)%ptr => this%scalars%scalar_fields(k)%s
       end if
    end do
    this%ev_cap%ptr => &
              neko_field_registry%get_field(fields(1+2*(1+this%n_scalars)))

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
       this%h_k%x(:,:,:,e) = diameter
    end do
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%h2%x, this%h2%x_d, this%u%dof%size(), &
                             HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%h_k%x, this%h_k%x_d, this%u%dof%size(), &
                             HOST_TO_DEVICE, sync = .false.)
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
     
     ! Compute the entropy at the first time step
     if (time%tstep .eq. 1) then
        call field_rzero(E_vel)
        if (this%if_filter) then
           ! filter the velocity magnitude all together
           call field_addcol3(E_vel, this%u, this%u)
           call field_addcol3(E_vel, this%v, this%v)
           call field_addcol3(E_vel, this%w, this%w)
           call field_sqrt(E_vel)
           call this%filter%apply(wa_vel, E_vel)
           call field_sub2(wa_vel, E_vel)
           call this%coef%gs_h%op(wa_vel, GS_OP_ADD)
           if (NEKO_BCKND_DEVICE .eq. 1) then
              call device_col2(wa_vel%x_d, coef%mult_d, n)
           else
              call col2(wa_vel%x, coef%mult, n)
           end if
           ! Here E_vel is just fu such that the residual is computed correctly
           call field_copy(E_vel, wa_vel)
        else
           call field_col3(wa_vel, this%u, this%u)
           call field_add2(E_vel, wa_vel)
           call field_col3(wa_vel, this%v, this%v)
           call field_add2(E_vel, wa_vel)
           call field_col3(wa_vel, this%w, this%w)
           call field_add2(E_vel, wa_vel)
           call field_sqrt(E_vel)
        end if
     end if

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
       ! Compute the entropy at the first time step
       if (time%tstep .eq. 1) then
          if (this%if_filter) then
             call this%filter%apply(wa_s_i, this%s(i)%ptr)
             call field_sub2(wa_s_i, this%s(i)%ptr)
             call this%coef%gs_h%op(wa_s_i, GS_OP_ADD)
             if (NEKO_BCKND_DEVICE .eq. 1) then
                call device_col2(wa_s_i%x_d, coef%mult_d, n)
             else
                call col2(wa_s_i%x, coef%mult, n)
             end if
             call field_copy(E_s_i, wa_s_i)
          else
             call field_copy(E_s_i, this%s(i)%ptr)
          end if
       end if
      
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
    type(field_t), pointer :: fu
    
    ! Arrays for scaling
    type(field_t), pointer :: E_true
    integer :: E_true_index

    type(field_t), pointer :: E_avg_field
    integer :: E_avg_index

    type(field_t), pointer :: E_var
    integer :: E_var_index

    type(field_t), pointer :: B_elem
    integer :: B_elem_index

    real(kind=rp) :: tol

    integer :: temp_index
    integer :: filt_field_index

    integer :: i, j, n
    real(kind=rp) :: scaling_factor
    
    ! create the work array and the upper bound
    call neko_scratch_registry%request_field(ta, temp_index, .false.)

    ! create the work array for filtered field
    call neko_scratch_registry%request_field(fu, filt_field_index, .false.)

    ! create work arrays for scaling
    call neko_scratch_registry%request_field(E_true, &
                                             E_true_index, .false.)

    call neko_scratch_registry%request_field(E_var, &
                                             E_var_index, .false.)

    call neko_scratch_registry%request_field(E_avg_field, &
                                             E_avg_index, .false.)

    call neko_scratch_registry%request_field(B_elem, B_elem_index, .false.)

    ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
    associate(u => this%u, v => this%v, w => this%w, E_vel => this%E(1), &
             ext_bdf => this%ext_bdf, &
             dt => time%dt, coef => this%coef, wa_vel => this%wa(1), &
             D_vel => this%D(1)%ptr, gs => this%coef%gs_h, &
             adv => this%adv, &
             Xh => this%coef%Xh, &
             ev_cap => this%ev_cap%ptr, &
             entropy_viscosity_vel => this%entropy_viscosity(1)%ptr)

    n = u%dof%size()

    if (this%if_filter) then
      ! filter the velocity magnitude all together
      call field_rzero(ta)
      call field_addcol3(ta, u, u)
      call field_addcol3(ta, v, v)
      call field_addcol3(ta, w, w)
      call field_sqrt(ta)
      call maxnorm_3d(ev_cap%x, ta%x, coef%Xh%lx, coef%msh%nelv)
      call this%filter%apply(fu, ta)
      call field_sub2(fu, ta)
      call gs%op(fu, GS_OP_ADD)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(fu%x_d, coef%mult_d, n)
      else
         call col2(fu%x, coef%mult, n)
      end if
      ! Here E_vel is just fu such that the residual is computed correctly
      call field_copy(E_vel, fu)
    else
      call field_col3(ta, u, u)
      call field_copy(E_vel, ta)
      call field_col3(ta, v, v)
      call field_add2(E_vel, ta)
      call field_col3(ta, w, w)
      call field_add2(E_vel, ta)
      call field_sqrt(E_vel)
      call maxnorm_3d(ev_cap%x, E_vel%x, coef%Xh%lx, coef%msh%nelv)
    end if

    call field_col2(ev_cap, this%h_k)
    ! now ev_cap is a work array for the upper bound of the viscosity
    call field_cmult(ev_cap, this%c_max)

    ! temporal derivative
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
   
    ! multiply 2 and the filtered field itself to get the real residual
    call field_col2(D_vel, fu)
    call field_cmult(D_vel, 2.0_rp)
    
    call field_copy(entropy_viscosity_vel, D_vel)
    call field_absval(entropy_viscosity_vel)

    ! Correct E_vel to be fu^2
    call field_col3(E_true, E_vel, E_vel)

    call dottnsr_3d(E_avg_field%x, E_true%x, &
         coef%B, coef%Xh%lx, coef%msh%nelv)
    call dot1tnsr_3d(B_elem%x, coef%B, coef%Xh%lx, coef%msh%nelv)
    call field_invcol2(E_avg_field, B_elem)
    call field_cmult(E_avg_field, -1.0_rp)

    call field_add3(E_var, E_true, E_avg_field)
    call field_absval(E_var)

    call field_cmult(entropy_viscosity_vel, &
         this%c_E)

    call maxnorm_3d(ta%x, E_var%x, coef%Xh%lx, coef%msh%nelv)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       tol = this%tol_coef * (device_glmax(E_true%x_d, n) - &
                              device_glmin(E_true%x_d, n))
    else
       tol = this%tol_coef * (glmax(E_true%x, n) - glmin(E_true%x, n))
    end if 
    call field_invcol2_nonzero(entropy_viscosity_vel, ta, tol)

    call field_col2(entropy_viscosity_vel, this%h2)
    call field_pwmin2(entropy_viscosity_vel, ev_cap)

    call gs%op(entropy_viscosity_vel, GS_OP_ADD)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(entropy_viscosity_vel%x_d, coef%mult_d, n)
    else
       call col2(entropy_viscosity_vel%x, coef%mult, n)
    end if

   end associate

    do i = 1, this%n_scalars
       ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
       associate(s_i => this%s(i)%ptr, E_s_i => this%E(i+1), &
                 ext_bdf => this%ext_bdf, &
                 dt => time%dt, coef => this%coef, wa_s_i => this%wa(i+1), &
                 D_s_i => this%D(i+1)%ptr, gs => this%coef%gs_h, &
                 adv => this%adv, &
                 u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
                 ev_cap => this%ev_cap%ptr, &
                 entropy_viscosity_i => this%entropy_viscosity(i+1)%ptr)

       n = s_i%dof%size()

       if (this%if_filter) then
         call this%filter%apply(fu, s_i)
         call field_sub2(fu, s_i)
         call gs%op(fu, GS_OP_ADD)
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_col2(fu%x_d, coef%mult_d, n)
         else
            call col2(fu%x, coef%mult, n)
         end if
         call field_copy(E_s_i, fu)
       else
         call field_copy(E_s_i, s_i)
       end if

       ! temporal derivative
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
       
       ! multiply 2 and the filtered field itself to get the real residual
       call field_col2(D_s_i, fu)
       call field_cmult(D_s_i, 2.0_rp)

       call field_copy(entropy_viscosity_i, D_s_i)
       call field_absval(entropy_viscosity_i)

       ! Correct E_vel to be fs^2
       call field_col3(E_true, E_s_i, E_s_i)

       call dottnsr_3d(E_avg_field%x, E_true%x, &
            coef%B, coef%Xh%lx, coef%msh%nelv)
       call dot1tnsr_3d(B_elem%x, coef%B, coef%Xh%lx, coef%msh%nelv)
       call field_invcol2(E_avg_field, B_elem)
       call field_cmult(E_avg_field, -1.0_rp)

       call field_add3(E_var, E_true, E_avg_field)
       call field_absval(E_var)

       call field_cmult(entropy_viscosity_i, &
            this%c_E)

       call maxnorm_3d(ta%x, E_var%x, coef%Xh%lx, coef%msh%nelv)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          tol = this%tol_coef * (device_glmax(E_true%x_d, n) - &
                                 device_glmin(E_true%x_d, n))
       else
          tol = this%tol_coef * (glmax(E_true%x, n) - glmin(E_true%x, n))
       end if 

       call field_invcol2_nonzero(entropy_viscosity_i, ta, tol)

       call field_col2(entropy_viscosity_i, this%h2)
       call field_pwmin2(entropy_viscosity_i, ev_cap)

       call gs%op(entropy_viscosity_i, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(entropy_viscosity_i%x_d, coef%mult_d, n)
       else
          call col2(entropy_viscosity_i%x, coef%mult, n)
       end if

       end associate
    end do

    call neko_scratch_registry%relinquish_field(temp_index)
    call neko_scratch_registry%relinquish_field(filt_field_index)
    call neko_scratch_registry%relinquish_field(E_true_index)
    call neko_scratch_registry%relinquish_field(E_var_index)
    call neko_scratch_registry%relinquish_field(E_avg_index)
    call neko_scratch_registry%relinquish_field(B_elem_index)

  end subroutine entropy_viscosity_compute

end module entropy_viscosity
