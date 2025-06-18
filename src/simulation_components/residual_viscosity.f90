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
!> A simulation component that computes residual_viscosity
!! The values are stored in the field registry under the name 'residual_viscosity'

module residual_viscosity
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use json_utils, only : json_get
  use field, only : field_t, field_ptr_t
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
  use field_math, only : field_col3, field_copy, field_absval, field_rzero, &
                         field_cmult, field_sub2, field_col2
  use math, only : invcol2, col2, cadd, absval, cfill, vlsc2
  use device_math, only : device_invcol2
  use gather_scatter, only : GS_OP_ADD
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: residual_viscosity_t
     !> coefficient
     real(kind=rp) :: c_E
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> Scalar field
     integer :: n_scalars = 0
     type(field_ptr_t), allocatable :: s(:)
     !> coef
     type(coef_t), pointer :: coef
     !> Some field to be used
     type(field_t), allocatable :: abx1(:), abx2(:)
     type(field_t) :: volume_el
     type(field_series_ptr_t), allocatable :: slag(:)

     !> X residual_viscosity component.
     type(field_ptr_t), allocatable :: residual_viscosity(:)

     !> Residual.
     type(field_t), allocatable :: D(:)
     !> work array.
     type(field_t), allocatable :: wa(:)

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
     procedure, pass(this) :: init => residual_viscosity_init_from_json
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => residual_viscosity_init_common
     !> Destructor.
     procedure, pass(this) :: free => residual_viscosity_free
     !> Part of the residual viscosity computation before the time stepping
     procedure, pass(this) :: preprocess_ => residual_viscosity_preprocess
     !> Part of the residual viscosity computation after the time stepping
     procedure, pass(this) :: compute_ => residual_viscosity_compute
  end type residual_viscosity_t

contains

  !> Constructor from json.
  subroutine residual_viscosity_init_from_json(this, json, case)
    class(residual_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case

    call this%init_base(json, case)
    
    call json_get(json, "c_E", this%c_E)

    call this%init_common(json, case)
  end subroutine residual_viscosity_init_from_json

  !> Common part of constructors.
  subroutine residual_viscosity_init_common(this, json, case)
    class(residual_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20), allocatable :: fields(:)
    integer :: e, k
    real(kind=rp) :: volume_element

    this%scalars => case%scalars
    this%n_scalars = size(this%scalars%scalar_fields)
    allocate(fields(this%n_scalars))
    do k = 1, this%n_scalars
       write(fields(k), '(A,I0)') 'res_visc_s', k
    end do
    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/residual_viscosity.f90\the registry.
    call json%add("fields", fields)
    call this%writer%init(json, case)

    this%coef => case%fluid%c_Xh

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
      this%makebdf => f1%makebdf
      this%makeext => f1%makeabf
      this%adv => f1%adv
      this%ext_bdf => f1%ext_bdf
    class default
      call neko_error("For fluid, residual &
      &viscosity currently only support pnpn scheme")
    end select

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    call this%volume_el%init(this%u%dof)

    allocate(this%slag(this%n_scalars))
    allocate(this%D(this%n_scalars))
    allocate(this%wa(this%n_scalars))
    allocate(this%abx1(this%n_scalars))
    allocate(this%abx2(this%n_scalars))
    allocate(this%s(this%n_scalars))
    allocate(this%residual_viscosity(this%n_scalars))

    do k = 1, this%n_scalars
       this%residual_viscosity(k)%ptr => &
              neko_field_registry%get_field(fields(k))

       call this%D(k)%init(this%u%dof)
       call this%wa(k)%init(this%u%dof)
       call this%abx1(k)%init(this%u%dof)
       call this%abx2(k)%init(this%u%dof)
       
       this%s(k)%ptr => this%scalars%scalar_fields(k)%s
       this%slag(k)%ptr => this%scalars%scalar_fields(k)%slag

       call field_rzero(this%wa(k))
    end do

    do e = 1, this%coef%msh%nelv
       volume_element = 0.0_rp
       do k = 1, this%coef%Xh%lx * this%coef%Xh%ly * this%coef%Xh%lz
          volume_element = volume_element + this%coef%B(k, 1, 1, e)
       end do
       this%volume_el%x(:,:,:,e) = volume_element
    end do

  end subroutine residual_viscosity_init_common

  !> Destructor.
  subroutine residual_viscosity_free(this)
    class(residual_viscosity_t), intent(inout) :: this
    call this%free_base()
  end subroutine residual_viscosity_free

  !> Part of the residual_viscosity computation before the time stepping.
  !! @param time The time state.
  subroutine residual_viscosity_preprocess(this, time)
    class(residual_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, n
    
    do i = 1, this%n_scalars
       !! estimate by the difference of the extrapolated and the solved advection
       associate(u => this%u, v => this%v, w => this%w, &
                 s => this%s(i)%ptr, D => this%D(i), &
                 abx1 => this%abx1(i), abx2 => this%abx2(i), &
                 coef => this%coef, &
                 rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
                 adv => this%adv, &
                 makeext => this%makeext, ext_bdf => this%ext_bdf, &
                 Xh => this%coef%Xh)
      
       n = s%dof%size()
       call field_rzero(D)
       call adv%compute_scalar(u, v, w, s, D, &
                Xh, coef, n)
       call makeext%compute_scalar(abx1, abx2, D%x, &
                rho%x(1,1,1,1), ext_bdf%advection_coeffs, n)

       end associate
       
      !  !! estimate by ds/dt + ui ds/dxi
      !  associate(s => this%s(i)%ptr, wa => this%wa(i), &
      !            coef => this%coef, &
      !            rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
      !            makebdf => this%makebdf, &
      !            slag => this%slag(i)%ptr, ext_bdf => this%ext_bdf)
      
      !  n = s%dof%size()
      !  call field_rzero(wa)
      !  call makebdf%compute_scalar(slag, wa%x, s, coef%B, rho%x(1,1,1,1), &
      !          dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
      !  end associate
    end do

  end subroutine residual_viscosity_preprocess

  !> Part of the residual_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine residual_viscosity_compute(this, time)
    class(residual_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ta ! temporal array
    integer :: temp_indices
    integer :: i, n

    do i = 1, this%n_scalars

       call neko_scratch_registry%request_field(ta, temp_indices)

       !! estimate by the difference of the extrapolated and the solved advection
       associate(s => this%s(i)%ptr, ext_bdf => this%ext_bdf, &
                 dt => time%dt, coef => this%coef, D => this%D(i), &
                 adv => this%adv, &
                 u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
                 residual_viscosity => this%residual_viscosity(i)%ptr)

       n = s%dof%size()
       call field_rzero(ta)
       call adv%compute_scalar(u, v, w, s, ta, Xh, coef, n)
       call field_sub2(D, ta, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(D%x_d, coef%B_d, n)
       else
          call invcol2(D%x, coef%B, n)
       end if
       call field_copy(residual_viscosity, D)
       call field_absval(residual_viscosity)
       ! it should be scaled by f(ext_bdf%diffusion_time_order) and also dt
       ! preliminary, f could be 0.01937*exp(-5.7363*ext_bdf%diffusion_time_order)
       ! Could be determined afterwards
       call field_cmult(residual_viscosity, &
           this%c_E)

       end associate

       !! estimate by ds/dt + ui ds/dxi
      !  associate(s => this%s(i)%ptr, ext_bdf => this%ext_bdf, &
      !            dt => time%dt, coef => this%coef, wa => this%wa(i), &
      !            D => this%D(i), gs => this%coef%gs_h, &
      !            adv => this%adv, &
      !            u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
      !            residual_viscosity => this%residual_viscosity(i)%ptr)

      !  n = s%dof%size()

      !  call field_copy(ta, s)
      !  call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
      !  if (NEKO_BCKND_DEVICE .eq. 1) then
      !     call device_invcol2(wa%x_d, coef%B_d, n)
      !  else
      !     call invcol2(wa%x, coef%B, n)
      !  end if
      !  call field_sub2(ta, wa)
      !  call field_copy(D, ta)

      !  ! advection part
      !  call field_rzero(ta, n)
      !  call adv%compute_scalar(u, v, w, s, ta, &
      !         Xh, coef, n)
      !  if (NEKO_BCKND_DEVICE .eq. 1) then
      !     call device_invcol2(ta%x_d, coef%B_d, n)
      !  else
      !     call invcol2(ta%x, coef%B, n)
      !  end if
      !  call field_sub2(D, ta, n)
      !  call field_copy(residual_viscosity, D)
      !  call field_absval(residual_viscosity)

      !  ! it should be scaled by f(ext_bdf%diffusion_time_order)
      !  ! preliminary, f could be 0.01937*exp(-5.7363*ext_bdf%diffusion_time_order)
      !  ! Could be determined afterwards
      !  call field_cmult(residual_viscosity, &
      !      this%c_E)

      !  end associate

       call neko_scratch_registry%relinquish_field(temp_indices)

    end do

  end subroutine residual_viscosity_compute

end module residual_viscosity
