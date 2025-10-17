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
  use field, only : field_t, field_ptr_t
  use utils, only : neko_error, neko_type_error
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use coefs, only : coef_t
  use math, only : cfill, copy, rzero
  use device_math, only : device_rzero, device_cfill, device_copy
  use device, only : device_map
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use neko_config, only : NEKO_BCKND_DEVICE
  implicit none
  private

  ! List of all possible directions for filtering
  character(len=20) :: KNOWN_DIRECTIONS(7) = [character(len=3) :: &
       "rst", &
       "rs", "rt", "st", &
       "r", "s", "t"]

  !> Implements the spectral vanishing viscosity.
  type, public :: svv_t
    !> filter
    type(elementwise_filter_t) :: filter
    !> filtering direction
    character(len=:), allocatable:: direction
    !> Power coefficient for the SVV kernel
    real(kind=rp) :: power_coef
    !> coef
    type(coef_t), pointer :: coef
    !> the viscosity field
    real(kind=rp), allocatable :: h1(:,:,:,:)
    type(c_ptr) :: h1_d = C_NULL_PTR
    !> a pointer pointing to a potentially variable viscosity field
    character(len=:), allocatable :: nue_field_name
    type(field_t), pointer :: nue
    type(field_ptr_t) :: nue_uvw(3)
    !> a logical to identify whether h1 is time variable
    logical :: tvar_h1 = .false.
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
    character(len=:), allocatable :: nu_type, nue_field_name_tmp, direction
    integer :: i

    this%coef => coef

    ! set up the viscosity coefficient field
    allocate(this%h1(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%h1, this%h1_d, this%coef%dof%size())
       call device_rzero(this%h1_d, this%coef%dof%size())
    else
       call rzero(this%h1, this%coef%dof%size())
    end if
    
    call json_get_or_default(json, "svv.direction", &
         direction, "rst")
    this%direction = trim(direction)
    if (this%direction .ne. "rst" .and. &
        this%direction .ne. "rs" .and. &
        this%direction .ne. "rt" .and. &
        this%direction .ne. "st" .and. &
        this%direction .ne. "r" .and. &
        this%direction .ne. "s" .and. &
        this%direction .ne. "t") then
       call neko_type_error("The direction of the SVV ", &
            this%direction, KNOWN_DIRECTIONS)
    end if

    call json_get(json, "svv.power_coefficient", this%power_coef)

    call json_get(json, "svv.nu.type", nu_type)
    select case (trim(nu_type))
    case ("value")
       call json_get(json, "svv.nu.value", nu_val)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_cfill(this%h1_d, nu_val, coef%dof%size())
       else
          call cfill(this%h1, nu_val, coef%dof%size())
       end if 
    case ("field")
       call json_get_or_default(json, "svv.nu.time_variable", this%tvar_h1, .true.)
       call json_get(json, "svv.nu.field_name", this%nue_field_name)
       if (neko_field_registry%field_exists(this%nue_field_name)) then
          this%nue => neko_field_registry%get_field(this%nue_field_name)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call device_copy(this%h1_d, this%nue%x_d, this%coef%dof%size())
          else
             call copy(this%h1, this%nue%x, this%coef%dof%size())
          end if
       end if
    case default
       call neko_error("Invalid nu.type for svv")
    end select

    ! set up the filter
    this%filter%filter_type = "nonBoyd"
    call this%filter%init_from_components(coef%Xh%lx)
    ! assign the SVV Kernel
    do i = 1, this%coef%Xh%lx
       this%filter%transfer(i) = ((i - 1.0_rp) / (this%coef%Xh%lx - 1.0_rp)) &
                              ** ((this%coef%Xh%lx - 1.0_rp) * this%power_coef)
       this%filter%transfer(i) = 1.0_rp - this%filter%transfer(i)
    end do
    ! build the 1d elementwise filter
    call this%filter%build_1d()

  end subroutine svv_init_from_json

  !> Update of h1 is it's time varying
  subroutine update_h1(this, tstep)
    class(svv_t), intent(inout) :: this
    integer, intent(in) :: tstep

    if (.not. this%tvar_h1) return

    if (tstep .eq. 1) then
       this%nue => neko_field_registry%get_field(this%nue_field_name)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%h1_d, this%nue%x_d, this%coef%dof%size())
    else
       call copy(this%h1, this%nue%x, this%coef%dof%size())
    end if

  end subroutine update_h1

end module spectral_vanishing_viscosity