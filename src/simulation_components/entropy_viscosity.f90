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
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use field, only : field_t
  use field_series, only : field_series_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use fluid_pnpn, only : fluid_pnpn_t
  use scalar_pnpn, only : scalar_pnpn_t
  use utils, only : neko_error
  use field_math, only : field_col3, field_copy, field_absval, field_rzero, &
                         field_cmult, field_sub2 
  use math, only : invcol2, col2
  use gather_scatter, only : GS_OP_ADD
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: entropy_viscosity_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w
     !> Scalar field
     type(field_t), pointer :: s
     !> Some field to be used, currently call it s2
     type(field_t) :: s2
     type(field_series_t) :: s2lag

     !> X entropy_viscosity component.
     type(field_t), pointer :: entropy_viscosity

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
     procedure, pass(this) :: init => entropy_viscosity_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          entropy_viscosity_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          entropy_viscosity_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => entropy_viscosity_init_common
     !> Destructor.
     procedure, pass(this) :: free => entropy_viscosity_free
     !> Part of the entropy viscosity computation before the time stepping
     procedure, pass(this) :: preprocess_ => entropy_viscosity_preprocess
     !> Part of the entropy viscosity computation after the time stepping
     procedure, pass(this) :: compute_ => entropy_viscosity_compute
  end type entropy_viscosity_t

contains

  !> Constructor from json.
  subroutine entropy_viscosity_init_from_json(this, json, case)
    class(entropy_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20) :: fields(1)
    type(field_t), pointer :: u, v, w, entropy_viscosity

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/entropy_viscosity.f90\the registry.
    fields(1) = "entropy_viscosity"
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call this%init_common(case)
  end subroutine entropy_viscosity_init_from_json

  !> Common part of constructors.
  subroutine entropy_viscosity_init_common(this, case)
    class(entropy_viscosity_t), intent(inout) :: this
    class(case_t), intent(inout), target ::case

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    this%s => neko_field_registry%get_field("s")
    this%entropy_viscosity => neko_field_registry%get_field("entropy_viscosity")

    call this%s2%init(this%u%dof)
    call this%s2lag%init(this%s2, 2)
    call this%D%init(this%u%dof)
    call this%wa%init(this%u%dof)

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
    class default
      call neko_error("For fluid, entropy &
      &viscosity currently only support pnpn scheme")
    end select

    this%scalar => case%scalar

  end subroutine entropy_viscosity_init_common

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine entropy_viscosity_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       filename, precision)
    class(entropy_viscosity_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)
    fields(1) = "entropy_viscosity"

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%writer%init_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller, fields, filename, precision)
    call this%init_common(case)

  end subroutine entropy_viscosity_init_from_controllers

  !> Constructor from components, passing properties to the
  !! time_based_controller` components in the base type.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine entropy_viscosity_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, filename, precision)
    class(entropy_viscosity_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)
    fields(1) = "entropy_viscosity"

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%writer%init_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, fields, filename, precision)
    call this%init_common(case)

  end subroutine entropy_viscosity_init_from_controllers_properties

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
    integer :: n

    associate(s => this%s, s2 => this%s2, wa => this%wa, &
              coef => this%scalar%c_Xh, &
              rho => this%scalar%rho, dt => time%dt, &
              makebdf => this%scalar%makebdf, &
              s2lag => this%s2lag, ext_bdf => this%fluid%ext_bdf)
    
    n = s%dof%size()

    call field_rzero(wa)
    call makebdf%compute_scalar(s2lag, wa%x, s2, coef%B, rho%x(1,1,1,1), &
            dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
    call s2lag%update()

    end associate

  end subroutine entropy_viscosity_preprocess

  !> Part of the entropy_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine entropy_viscosity_compute(this, time)
    class(entropy_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ta ! temporal array
    integer :: temp_indices(1)
    integer :: n

    call neko_scratch_registry%request_field(ta, temp_indices(1))

    associate(s => this%s, s2 => this%s2, ext_bdf => this%fluid%ext_bdf, &
              dt => time%dt, coef => this%scalar%c_Xh, wa => this%wa, &
              D => this%D, gs => this%scalar%gs_Xh, adv => this%scalar%adv, &
              u => this%u, v => this%v, w => this%w, Xh => this%scalar%Xh, &
              entropy_viscosity => this%entropy_viscosity)

    n = s%dof%size()
    call field_rzero(s2, n)

    ! local temporal change part
    ! call field_copy(entropy_viscosity, wa)
    !  call field_col3(s2, s, s, n)
    call field_copy(s2, s, n)
    !  call field_absval(s2)
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
    call gs%op(ta, GS_OP_ADD)
    call col2(ta%x, coef%mult, n)
    call field_sub2(D, ta, n)
    call field_copy(entropy_viscosity, D)
    call field_absval(entropy_viscosity)
    ! it should be scaled by f(ext_bdf%diffusion_time_order)
    ! preliminary, f = 0.01937*exp(-5.7363*ext_bdf%diffusion_time_order)
    ! Could be determined afterwards
    call field_cmult(entropy_viscosity, &
         1.0_rp/dt**ext_bdf%diffusion_time_order)
    ! call field_col2(entropy_viscosity, s, n)
    ! call field_absval(entropy_viscosity)

    end associate

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine entropy_viscosity_compute

end module entropy_viscosity
