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
  use json_utils, only : json_get, json_extract_object
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
                         field_invcol2, field_sqrt, field_add2
  use math, only : invcol2, col2, glsum, glsc2, glmax
  use device_math, only : device_invcol2, device_col2, device_glsum, &
                          device_glsc2, device_glmax
  use gather_scatter, only : GS_OP_ADD
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: entropy_viscosity_t
     !> coefficient
     real(kind=rp) :: c_E
     !> A low pass filter for the field
     type(elementwise_filter_t) :: filter
     logical :: if_filter = .false.
     !> X velocity component.
     type(field_t), pointer :: u
     type(field_t) :: E_vel_var
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> Velocity magnitude.
     type(field_t) :: vel_mag
     !> work array for the temporal derivative
     type(field_t), allocatable :: wa(:)
     !> Scalar field
     integer :: n_scalars = 0
     type(field_ptr_t), allocatable :: s(:)
     type(field_t), allocatable :: E_s_var(:)
     type(field_t), allocatable :: E(:)
     type(field_series_t), allocatable :: Elag(:)
     !> coef
     type(coef_t), pointer :: coef
     !> Some field to be used
     type(field_t), allocatable :: abx1(:), abx2(:)
     type(field_t) :: h2
     real(kind=rp) :: volume_domain

     !> X entropy_viscosity component.
     type(field_ptr_t), allocatable :: entropy_viscosity(:)

     !> Residual.
     type(field_t), allocatable :: D(:)

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

    call this%init_base(json, case)
    
    call json_get(json, "c_E", this%c_E)

    call this%init_common(json, case)
  end subroutine entropy_viscosity_init_from_json

  !> Common part of constructors.
  subroutine entropy_viscosity_init_common(this, json, case)
    class(entropy_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20), allocatable :: fields(:)
    integer :: e, k
    real(kind=rp) :: volume_element

    if (allocated(case%scalars)) then
       this%scalars => case%scalars
       this%n_scalars = size(this%scalars%scalar_fields)
    else
       this%n_scalars = 0
    end if
    allocate(fields(1+this%n_scalars))
    fields(1) = 'entr_visc_vel'
    do k = 1, this%n_scalars
       write(fields(k+1), '(A,I0)') 'entr_visc_s', k
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

    call this%E_vel_var%init(this%u%dof)
    call this%vel_mag%init(this%u%dof)

    call this%h2%init(this%u%dof)

    if (this%n_scalars .ne. 0) then
       allocate(this%s(this%n_scalars))
       allocate(this%E_s_var(this%n_scalars))
    end if
    
    allocate(this%E(1+this%n_scalars))
    allocate(this%Elag(1+this%n_scalars))
    allocate(this%wa(1+this%n_scalars))

    allocate(this%D(1+this%n_scalars))
    allocate(this%abx1(1+this%n_scalars))
    allocate(this%abx2(1+this%n_scalars))
    
    allocate(this%entropy_viscosity(1+this%n_scalars))

    do k = 1, 1+this%n_scalars
       this%entropy_viscosity(k)%ptr => &
              neko_field_registry%get_field(fields(k))

       call this%E(k)%init(this%u%dof)
       call this%Elag(k)%init(this%E(k), 2)
       call this%wa(k)%init(this%u%dof)

       call this%D(k)%init(this%u%dof)
       call this%abx1(k)%init(this%u%dof)
       call this%abx2(k)%init(this%u%dof)
       
       if (k .le. this%n_scalars) then
          this%s(k)%ptr => this%scalars%scalar_fields(k)%s
          call this%E_s_var(k)%init(this%u%dof)
       end if
    end do

    do e = 1, this%coef%msh%nelv
       volume_element = 0.0_rp
       do k = 1, this%coef%Xh%lx * this%coef%Xh%ly * this%coef%Xh%lz
          volume_element = volume_element + this%coef%B(k, 1, 1, e)
       end do
       this%h2%x(:,:,:,e) = volume_element**(1.0_rp/3.0_rp) * &
                            volume_element**(1.0_rp/3.0_rp) / &
                            (this%coef%Xh%lx-1.0_rp) / &
                            (this%coef%Xh%lx-1.0_rp)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%h2%x, this%h2%x_d, this%u%dof%size(), &
                          HOST_TO_DEVICE, sync = .false.)
       this%volume_domain = device_glsum(this%coef%B_d, this%u%dof%size())
    else
       this%volume_domain = glsum(this%coef%B, this%u%dof%size())
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
    associate(wa => this%wa(1), E => this%E(1), Elag => this%Elag(1), &
              coef => this%coef, &
              rho => this%fluid%rho, dt => time%dt, &
              makebdf => this%makebdf, &
              ext_bdf => this%ext_bdf)

     n = wa%dof%size()
     call field_rzero(wa)
     call makebdf%compute_scalar(Elag, wa%x, E, coef%B, &
               rho%x(1,1,1,1), dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
     call Elag%update()
     
    end associate
    
    do i = 1, this%n_scalars
       ! Time lag part for the BDF scheme of dE/dt
       associate(wa => this%wa(i+1), &
                 coef => this%coef, &
                 rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
                 makebdf => this%makebdf, E => this%E(i+1), &
                 Elag => this%Elag(i+1), ext_bdf => this%ext_bdf)
      
       call field_rzero(wa)
       call makebdf%compute_scalar(Elag, wa%x, E, coef%B, rho%x(1,1,1,1), &
               dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(wa%x_d, coef%B_d, n)
       else
          call invcol2(wa%x, coef%B, n)
       end if
       call Elag%update()

       end associate
    end do

  end subroutine entropy_viscosity_preprocess

  !> Part of the entropy_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine entropy_viscosity_compute(this, time)
    class(entropy_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_ptr_t) :: ta(1+this%n_scalars) ! temporal array
    type(field_ptr_t) :: fs(this%n_scalars)
    type(field_t), pointer :: fu, fv, fw
    real(kind=rp) :: E_vel_avg
    real(kind=rp) :: E_s_avg(this%n_scalars)
    integer :: temp_indices(1+this%n_scalars)
    integer :: filt_field_indices(3+this%n_scalars)
    integer :: i, j, n
    real(kind=rp) :: scaling_vel, scaling_s(this%n_scalars)

    do i = 1, 1+this%n_scalars
       call neko_scratch_registry%request_field(ta(i)%ptr, temp_indices(i))
    end do

    call neko_scratch_registry%request_field(fu, filt_field_indices(1))
    call neko_scratch_registry%request_field(fv, filt_field_indices(2))
    call neko_scratch_registry%request_field(fw, filt_field_indices(3))
    do i = 1, this%n_scalars
       call neko_scratch_registry%request_field(fs(i)%ptr, &
            filt_field_indices(3+i))
    end do

    ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
    associate(u => this%u, v => this%v, w => this%w, E => this%E(1), ta => ta(1)%ptr, &
             ext_bdf => this%ext_bdf, &
             dt => time%dt, coef => this%coef, wa => this%wa(1), &
             D => this%D(1), gs => this%coef%gs_h, &
             adv => this%adv, &
             Xh => this%coef%Xh, &
             entropy_viscosity => this%entropy_viscosity(1)%ptr)

    n = u%dof%size()

    if (this%if_filter) then

      call this%filter%apply(fu, u)
      call this%filter%apply(fv, v)
      call this%filter%apply(fw, w)
      
      call field_sub2(fu, u)
      call field_sub2(fv, v)
      call field_sub2(fw, w)

      call gs%op(fu, GS_OP_ADD)
      call gs%op(fv, GS_OP_ADD)
      call gs%op(fw, GS_OP_ADD)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_col2(fu%x_d, coef%mult_d, n)
         call device_col2(fv%x_d, coef%mult_d, n)
         call device_col2(fw%x_d, coef%mult_d, n)
      else
         call col2(fu%x, coef%mult, n)
         call col2(fv%x, coef%mult, n)
         call col2(fw%x, coef%mult, n)
      end if

      call field_col3(ta, fu, fu)
      call field_copy(E, ta)
      call field_col3(ta, fv, fv)
      call field_add2(E, ta)
      call field_col3(ta, fw, fw)
      call field_add2(E, ta)
      call field_sqrt(E)
      
    else
      call field_col3(ta, u, u)
      call field_copy(E, ta)
      call field_col3(ta, v, v)
      call field_add2(E, ta)
      call field_col3(ta, w, w)
      call field_add2(E, ta)
      call field_sqrt(E)
    end if

    call field_copy(ta, E)
    call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
    call field_sub2(ta, wa)
    call field_copy(D, ta)

    ! advection part
    call field_rzero(ta, n)
    call adv%compute_scalar(u, v, w, E, ta, &
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
    call field_sub2(D, ta, n)
    call field_copy(entropy_viscosity, D)
    call field_absval(entropy_viscosity)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       E_vel_avg = - device_glsc2(E%x_d, coef%B_d, n) / this%volume_domain
    else
       E_vel_avg = - glsc2(E%x, coef%B, n) / this%volume_domain
    end if
    call field_cadd2(this%E_vel_var, E, E_vel_avg)
    call field_absval(this%E_vel_var)
   
    if (NEKO_BCKND_DEVICE .eq. 1) then
       scaling_vel = this%c_E / device_glmax(this%E_vel_var%x_d, u%dof%size())
    else
       scaling_vel = this%c_E / glmax(this%E_vel_var%x, u%dof%size())
    end if

    call field_cmult(entropy_viscosity, &
         scaling_vel)
    call field_col2(entropy_viscosity, this%h2)

   end associate

    do i = 1, this%n_scalars
       ! The updated part for the BDF scheme of dE/dt and the updated ui dE/dxi
       associate(s => this%s(i)%ptr, &
                 fs => fs(i)%ptr, E => this%E(i+1), ta => ta(i+1)%ptr, &
                 ext_bdf => this%ext_bdf, &
                 dt => time%dt, coef => this%coef, wa => this%wa(i+1), &
                 D => this%D(i+1), gs => this%coef%gs_h, &
                 adv => this%adv, &
                 u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
                 E_s_avg => E_s_avg(i), E_s_var => this%E_s_var(i), &
                 scaling_s => scaling_s(i), &
                 entropy_viscosity => this%entropy_viscosity(i+1)%ptr)

       n = s%dof%size()

       if (this%if_filter) then
         call this%filter%apply(fs, s)
         call field_sub2(fs, s)
         call gs%op(fs, GS_OP_ADD)
         if (NEKO_BCKND_DEVICE .eq. 1) then
            call device_col2(fs%x_d, coef%mult_d, n)
         else
            call col2(fs%x, coef%mult, n)
         end if
         call field_copy(E, fs)
       else
         call field_copy(E, s)
       end if

       call field_absval(E)
       call field_copy(ta, E)
       call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
       call field_sub2(ta, wa)
       call field_copy(D, ta)

       ! advection part
       call field_rzero(ta, n)
       call adv%compute_scalar(u, v, w, E, ta, &
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
       call field_sub2(D, ta, n)
       call field_copy(entropy_viscosity, D)
       call field_absval(entropy_viscosity)
       
       if (NEKO_BCKND_DEVICE .eq. 1) then
          E_s_avg = - device_glsc2(E%x_d, coef%B_d, n) / this%volume_domain
       else
          E_s_avg = - glsc2(E%x, coef%B, n) / this%volume_domain
       end if
       call field_cadd2(E_s_var, E, E_s_avg)
       call field_absval(E_s_var)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          scaling_s = this%c_E / device_glmax(E_s_var%x_d, u%dof%size())
       else
          scaling_s = this%c_E / glmax(E_s_var%x, u%dof%size())
       end if

       call field_cmult(entropy_viscosity, &
            scaling_s)
       call field_col2(entropy_viscosity, this%h2)

       end associate
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
    call neko_scratch_registry%relinquish_field(filt_field_indices)

  end subroutine entropy_viscosity_compute

end module entropy_viscosity
