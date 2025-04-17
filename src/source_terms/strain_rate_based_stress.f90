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
!> Implements the `strain_rate_based_stress_t` type.
!! Maintainer: Shiyu Du.

module strain_rate_based_stress
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use strain_rate_based_stress_cpu, only : strain_rate_based_stress_compute_cpu
  implicit none
  private

  !> This source term adds the Coriolis force.
  type, public, extends(source_term_t) :: strain_rate_based_stress_t
     !> The "viscosity" term
     type(field_t), pointer :: nut
     !> The name for the nut field
     character(len=:), allocatable :: nut_field_name
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => strain_rate_based_stress_init_from_json
     !> The costrucructor from type components.
     procedure, pass(this) :: init_from_compenents => &
       strain_rate_based_stress_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => strain_rate_based_stress_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => strain_rate_based_stress_compute
  end type strain_rate_based_stress_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine strain_rate_based_stress_init_from_json(this, json, fields, coef)
    class(strain_rate_based_stress_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp) :: start_time, end_time
    character(len=:), allocatable :: nut_field_name

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    if (json%valid_path("nut_field")) then
       call json_get(json, "nut_field", nut_field_name)
       this%nut_field_name = nut_field_name
    else
       call neko_error("Specify the name for the nut_field")
    end if

    call strain_rate_based_stress_init_from_components(this, fields, &
         coef, start_time, end_time)

  end subroutine strain_rate_based_stress_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param omega The rotation vector.
  !! @param u_geo The geostrophic wind.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine strain_rate_based_stress_init_from_components(this, fields, &
               coef, start_time, end_time)
    class(strain_rate_based_stress_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (fields%size() .ne. 3) then
       call neko_error("Number of fields for the strain rate based stresses &
       &must be 3.")
    end if
  end subroutine strain_rate_based_stress_init_from_components

  !> Destructor.
  subroutine strain_rate_based_stress_free(this)
    class(strain_rate_based_stress_t), intent(inout) :: this

    call this%free_base()
  end subroutine strain_rate_based_stress_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine strain_rate_based_stress_compute(this, t, tstep)
    class(strain_rate_based_stress_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
  
    this%nut => neko_field_registry%get_field(this%nut_field_name)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("The strain rate based stress &
       &is only implemented on the CPU")
    else
       call strain_rate_based_stress_compute_cpu(this%fields, this%nut, &
                                                 this%coef)
    end if
  end subroutine strain_rate_based_stress_compute

end module strain_rate_based_stress
