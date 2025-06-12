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
!> Implements `svv_t`.
module spectral_vanishing_viscosity
  use num_types, only : rp
  use elementwise_filter, only: elementwise_filter_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use utils, only : neko_error
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use coefs, only : coef_t
  use math, only : cfill, copy
  implicit none
  private

  !> Implements the spectral vanishing viscosity.
  type, public :: svv_t
    !> filter
    type(elementwise_filter_t) :: filter
    !> coef
    type(coef_t), pointer :: coef
    !> the viscosity field
    real(kind=rp), allocatable :: h1(:,:,:,:)
    !> a pointer pointing to a potentially variable viscosity field
    type(field_t), pointer :: nue
    !> a logical to identify whether h1 is time variable
    logical :: tvar_h1
  contains
    procedure, pass(this) :: init => svv_init_from_json
    ! procedure, pass(this) :: free => svv_free
    procedure, pass(this) :: update => update_h1
  end type svv_t

contains
  !> Constructor
  subroutine svv_init_from_json(this, json, coef)
    class(svv_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in), target :: coef
    real(kind=rp) :: nu_val
    character(len=:), allocatable :: nu_type, nue_field_name
    integer :: i

    this%coef => coef

    ! set up the viscosity coefficient field
    allocate(this%h1(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    call json_get(json, "svv.nu.type", nu_type)
    select case (trim(nu_type))
    case ("value")
       call json_get(json, "svv.nu.value", nu_val)
       call cfill(this%h1, nu_val, coef%dof%size())
    case ("field")
       call json_get_or_default(json, "svv.nu.time_variable", this%tvar_h1, .true.)
       call json_get(json, "svv.nu.field_name", nue_field_name)
       this%nue => neko_field_registry%get_field(nue_field_name)
       call copy(this%h1, this%nue%x, coef%dof%size())
    case default
       call neko_error("Invalid nu.type for svv")
    end select

    ! set up the filter
    this%filter%filter_type = "nonBoyd"
    call this%filter%init_from_components(coef%Xh%lx, this%filter%filter_type)
    ! assign the SVV Kernel
    do i = 1, this%coef%Xh%lx
       this%filter%trnsfr(i) = ((i - 1.0_rp) / (this%coef%Xh%lx - 1.0_rp)) &
                              ** ((this%coef%Xh%lx - 1.0_rp) / 2.0_rp)
    end do
    ! build the 1d elementwise filter
    call this%filter%build_1d()

  end subroutine svv_init_from_json

  !> Update of h1 is it's time varying
  subroutine update_h1(this)
    class(svv_t), intent(inout) :: this

    call copy(this%h1, this%nue%x, this%coef%dof%size())

  end subroutine update_h1

end module spectral_vanishing_viscosity