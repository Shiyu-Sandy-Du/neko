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
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use json_utils, only : json_get
  use field, only : field_t
  use coefs, only : coef_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use fluid_pnpn, only : fluid_pnpn_t
  use scalar_pnpn, only : scalar_pnpn_t
  use utils, only : neko_error
  use field_math, only : field_col3, field_copy, field_absval, field_rzero, &
                         field_cmult, field_sub2, field_col2
  use math, only : invcol2, col2, cadd, absval, cfill, vlsc2
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
     type(field_t), pointer :: s
     !> coef
     type(coef_t), pointer :: coef
     !> Some field to be used
     type(field_t) :: abx1, abx2, volume_el
     type(field_t) :: s2, abs_var_s2
     type(field_series_t) :: s2lag

     !> X residual_viscosity component.
     type(field_t), pointer :: residual_viscosity

     !> Residual.
     type(field_t) :: D
     !> work array.
     type(field_t) :: wa

     !> Output writer.
     type(field_writer_t) :: writer

     !> A pointer point to the time scheme
     type(fluid_pnpn_t), pointer :: fluid
     type(scalar_pnpn_t), pointer :: scalar

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
    character(len=20) :: fields(1)
    type(field_t), pointer :: u, v, w, residual_viscosity

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/residual_viscosity.f90\the registry.
    fields(1) = "residual_viscosity"
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call json_get(json, "c_E", this%c_E)

    call this%init_common(case)
  end subroutine residual_viscosity_init_from_json

  !> Common part of constructors.
  subroutine residual_viscosity_init_common(this, case)
    class(residual_viscosity_t), intent(inout) :: this
    class(case_t), intent(inout), target ::case
    integer :: e, k
    real(kind=rp) :: volume_element

    this%coef => case%fluid%c_Xh

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    this%s => neko_field_registry%get_field("s")
    this%residual_viscosity => neko_field_registry%get_field("residual_viscosity")

    call this%s2%init(this%u%dof)
    call this%abs_var_s2%init(this%u%dof)
    call this%s2lag%init(this%s2, 2)
    call this%D%init(this%u%dof)
    call this%wa%init(this%u%dof)
    call this%abx1%init(this%u%dof)
    call this%abx2%init(this%u%dof)
    call this%volume_el%init(this%u%dof)

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
    class default
      call neko_error("For fluid, residual &
      &viscosity currently only support pnpn scheme")
    end select

    this%scalar => case%scalar
    call field_rzero(this%wa)

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
    integer :: n

    !!! estimate by the difference of the extrapolated and the solved advection
    ! associate(u => this%u, v => this%v, w => this%w, &
    !           s => this%s, wa => this%wa, &
    !           abx1 => this%abx1, abx2 => this%abx2, &
    !           coef => this%coef, &
    !           rho => this%scalar%rho, dt => time%dt, &
    !           adv => this%scalar%adv, &
    !           makeext => this%scalar%makeext, &
    !           s2lag => this%s2lag, ext_bdf => this%fluid%ext_bdf, &
    !           Xh => this%scalar%Xh)
    
    ! n = s%dof%size()
    ! call field_rzero(wa)
    ! call adv%compute_scalar(u, v, w, s, wa, &
    !           Xh, coef, n)
    ! call makeext%compute_scalar(abx1, abx2, wa%x, &
    !           rho%x(1,1,1,1), es2 => this%s2, xt_bdf%advection_coeffs, n)

    ! end associate
    
    !! estimate by ds/dt + ui ds/dxi
    associate(s => this%s, s2 => this%s2, wa => this%wa, &
              coef => this%coef, &
              rho => this%scalar%rho, dt => time%dt, &
              makebdf => this%scalar%makebdf, &
              s2lag => this%s2lag, ext_bdf => this%fluid%ext_bdf)
    
    n = s%dof%size()
    call field_rzero(wa)
    call makebdf%compute_scalar(s2lag, wa%x, s2, coef%B, rho%x(1,1,1,1), &
            dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
    call s2lag%update()

    end associate

  end subroutine residual_viscosity_preprocess

  !> Part of the residual_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine residual_viscosity_compute(this, time)
    class(residual_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ta ! temporal array
    integer :: temp_indices
    integer :: n, n_el, e
    real(kind=rp) :: int_s2_el, avg_s2_el
    real(kind=rp) :: tol = 1e-7

    ! !!! estimate by the difference of the extrapolated and the solved advection
    ! associate(u => this%u, v => this%v, w => this%w, &
    !           s => this%s, wa => this%wa, &
    !           D => this%D, &
    !           coef => this%coef, &
    !           rho => this%scalar%rho, dt => time%dt, &
    !           adv => this%scalar%adv, &
    !           s2lag => this%s2lag, &
    !           residual_viscosity => this%residual_viscosity, &
    !           ext_bdf => this%fluid%ext_bdf, Xh => this%scalar%Xh, &
    !           gs => this%scalar%gs_Xh)
    
    ! n = s%dof%size()
    ! call field_rzero(D)
    ! call adv%compute_scalar(u, v, w, s, D, &
    !           Xh, coef, n)
    ! call field_sub2(D, wa, n)
    ! call invcol2(D%x, coef%B, n)
    ! call gs%op(D, GS_OP_ADD)
    ! call col2(D%x, coef%mult, n)
    ! call field_copy(residual_viscosity, D)
    ! call field_absval(residual_viscosity)
    ! call field_cmult(residual_viscosity, &
    !      1.0_rp/dt**ext_bdf%diffusion_time_order)

    ! end associate

    !! estimate by ds/dt + ui ds/dxi
    call neko_scratch_registry%request_field(ta, temp_indices)

    associate(s => this%s, s2 => this%s2, ext_bdf => this%fluid%ext_bdf, &
              dt => time%dt, coef => this%coef, wa => this%wa, &
              D => this%D, gs => this%coef%gs_h, adv => this%scalar%adv, &
              u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
              residual_viscosity => this%residual_viscosity, &
              abs_var_s2 => this%abs_var_s2)

    n = s%dof%size()
    n_el = coef%Xh%lx*coef%Xh%ly*coef%Xh%lz

    ! call field_rzero(s2, n)

    ! local temporal change part
    call field_copy(s2, s, n)
    ! call field_absval(s2)
    ! call field_col3(s2, s, s, n)
    call field_copy(ta, s2)
    call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
    call invcol2(wa%x, coef%B, n)
    call field_sub2(ta, wa)
    call field_copy(D, ta)

    ! advection part
    call field_rzero(ta, n)
    call adv%compute_scalar(u, v, w, s2, ta, &
            Xh, coef, n)
    call invcol2(ta%x, coef%B, n)
    ! call gs%op(ta, GS_OP_ADD)
    ! call col2(ta%x, coef%mult, n)
    call field_sub2(D, ta, n)
    call field_copy(residual_viscosity, D)
    call field_absval(residual_viscosity)

    call field_copy(abs_var_s2, s2)
    do e = 1, coef%msh%nelv
       int_s2_el = vlsc2(s2%x(:,:,:,e), coef%B(:,:,:,e), n_el)
       avg_s2_el = -1.0_rp * int_s2_el/this%volume_el%x(1,1,1,e)
       call cadd(abs_var_s2%x(:,:,:,e), avg_s2_el, n_el)
       call absval(abs_var_s2%x(:,:,:,e), n_el)
       if (maxval(abs_var_s2%x(:,:,:,e)) .lt. tol) then
         call cfill(abs_var_s2%x(:,:,:,e), 0.0_rp, n_el)
       else
         call cfill(abs_var_s2%x(:,:,:,e), &
              1.0/maxval(abs_var_s2%x(:,:,:,e)), n_el)
       end if
    end do
    ! call gs%op(abs_var_s2, GS_OP_ADD)
    ! call col2(abs_var_s2%x, coef%mult, n)

    ! it should be scaled by f(ext_bdf%diffusion_time_order)
    ! preliminary, f = 0.01937*exp(-5.7363*ext_bdf%diffusion_time_order)
    ! Could be determined afterwards
    call field_cmult(residual_viscosity, &
         this%c_E)
    ! call field_col2(residual_viscosity, abs_var_s2)
    

    end associate

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine residual_viscosity_compute

end module residual_viscosity
