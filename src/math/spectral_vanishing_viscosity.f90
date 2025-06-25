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
  use utils, only : neko_error
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use coefs, only : coef_t
  use math, only : cfill, copy, rzero
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
    real(kind=rp), allocatable :: h1_1(:,:,:,:)
    real(kind=rp), allocatable :: h1_2(:,:,:,:)
    real(kind=rp), allocatable :: h1_3(:,:,:,:)
    !> a pointer pointing to a potentially variable viscosity field
    character(len=:), allocatable :: nue_field_name_1
    character(len=:), allocatable :: nue_field_name_2
    character(len=:), allocatable :: nue_field_name_3
    type(field_t), pointer :: nue
    type(field_ptr_t) :: nue_uvw(3)
    !> a logical to identify whether h1 is time variable
    logical :: tvar_h1 = .false.
    !> number of equation coupled
    integer :: eqn_number
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
    character(len=:), allocatable :: nu_type, nue_field_name_tmp
    integer :: i

    this%coef => coef

    ! set up the viscosity coefficient field
    allocate(this%h1(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(this%h1_1(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(this%h1_2(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    allocate(this%h1_3(coef%Xh%lx, coef%Xh%lx, coef%Xh%lx, coef%msh%nelv))
    call rzero(this%h1, this%coef%dof%size())
    call rzero(this%h1_1, this%coef%dof%size())
    call rzero(this%h1_2, this%coef%dof%size())
    call rzero(this%h1_3, this%coef%dof%size())
    call json_get(json, "svv.nu.type", nu_type)
    select case (trim(nu_type))
    case ("value")
       call json_get(json, "svv.nu.value", nu_val)
       call cfill(this%h1, nu_val, coef%dof%size())
       this%eqn_number = 1
    case ("field")
       call json_get_or_default(json, "svv.nu.time_variable", this%tvar_h1, .true.)
       call json_get(json, "svv.eqn_number", this%eqn_number)

       select case (this%eqn_number)
       case(1)
          call json_get(json, "svv.nu.field_name", this%nue_field_name_1)
       case(3)
          call json_get(json, "svv.nu.field_name", nue_field_name_tmp)
          allocate(character(len=len_trim(nue_field_name_tmp) + 2) &
                  :: this%nue_field_name_1)
          allocate(character(len=len_trim(nue_field_name_tmp) + 2) &
                  :: this%nue_field_name_2)
          allocate(character(len=len_trim(nue_field_name_tmp) + 2) &
                  :: this%nue_field_name_3)
          write(this%nue_field_name_1, '(A,A)') nue_field_name_tmp, "_u"
          write(this%nue_field_name_2, '(A,A)') nue_field_name_tmp, "_v"
          write(this%nue_field_name_3, '(A,A)') nue_field_name_tmp, "_w"
       case default
          call neko_error("SVV eqn_number has to be 1 or 3")
       end select

       
    case default
       call neko_error("Invalid nu.type for svv")
    end select

    ! set up the filter
    this%filter%filter_type = "nonBoyd"
    call this%filter%init_from_components(coef%Xh%lx, this%filter%filter_type)
    ! assign the SVV Kernel
    do i = 1, this%coef%Xh%lx
       this%filter%trnsfr(i) = ((i - 1.0_rp) / (this%coef%Xh%lx - 1.0_rp)) &
                              ** ((this%coef%Xh%lx - 1.0_rp) / 10.0_rp)
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
       if (this%eqn_number .eq. 1) then
          this%nue => neko_field_registry%get_field(this%nue_field_name_1)
       else
          this%nue_uvw(1)%ptr => &
            neko_field_registry%get_field(this%nue_field_name_1)
          this%nue_uvw(2)%ptr => &
            neko_field_registry%get_field(this%nue_field_name_2)
          this%nue_uvw(3)%ptr => &
            neko_field_registry%get_field(this%nue_field_name_3)
       end if
    end if

    if (this%eqn_number .eq. 1) then
       call copy(this%h1, this%nue%x, this%coef%dof%size())
    else
       call copy(this%h1_1, this%nue_uvw(1)%ptr%x, this%coef%dof%size())
       call copy(this%h1_2, this%nue_uvw(2)%ptr%x, this%coef%dof%size())
       call copy(this%h1_3, this%nue_uvw(3)%ptr%x, this%coef%dof%size())
    end if

  end subroutine update_h1

end module spectral_vanishing_viscosity