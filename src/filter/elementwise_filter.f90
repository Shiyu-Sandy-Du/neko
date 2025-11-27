! Copyright (c) 2024, The Neko Authors
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
!> Implements `elementwise_filter_t`.
module elementwise_filter
  use num_types, only : rp
  use filter, only: filter_t
  use math, only : rzero, rone, copy
  use field, only : field_t
  use coefs, only : coef_t
  use utils, only : neko_error, neko_type_error
  use neko_config, only : NEKO_BCKND_DEVICE
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get
  use speclib, only : zwgll, legendre_poly
  use matrix, only : matrix_t
  use mxm_wrapper, only : mxm
  use tensor, only : tnsr3d, trsp
  use device, only : device_map, device_free, device_memcpy, HOST_TO_DEVICE
  use device_math, only : device_cfill, device_rzero, device_glmax
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR, c_associated
  implicit none
  private
  
  ! List of all possible directions for filtering
  character(len=20) :: KNOWN_DIRECTIONS(7) = [character(len=3) :: &
       "rst", &
       "rs", "rt", "st", &
       "r", "s", "t"]

  !> Implements the elementwise filter for SEM.
  type, public, extends(filter_t) :: elementwise_filter_t
     !> filter type:
     !> possible options: "Boyd", "nonBoyd"
     character(len=:), allocatable :: filter_type
     !> dimension
     integer :: nx
     !> filtered wavenumber
     integer :: nt
     !> filtering direction
     character(len=:), allocatable :: direction
     !> matrix for 1d elementwise filtering
     real(kind=rp), allocatable :: fh(:,:), fht(:,:), ident(:,:)
     type(c_ptr) :: fh_d = C_NULL_PTR
     type(c_ptr) :: fht_d = C_NULL_PTR
     type(c_ptr) :: ident_d = C_NULL_PTR
     !> transfer function
     real(kind=rp), allocatable :: transfer(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => elementwise_filter_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_components => &
          elementwise_filter_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => elementwise_filter_free
     !> Set up 1D filter inside an element.
     procedure, pass(this) :: build_1d
     !> Filter a 3D field
     procedure, pass(this) :: apply => elementwise_field_filter_3d
  end type elementwise_filter_t

contains
  !> Constructor
  subroutine elementwise_filter_init_from_json(this, json, coef)
    class(elementwise_filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    real(kind=rp), allocatable :: transfer(:)
    character(len=:), allocatable :: filter_type
    character(len=:), allocatable :: direction

    ! Filter assumes lx = ly = lz
    call this%init_base(json, coef)

    call this%init_from_components(coef%dof%xh%lx)

    call json_get_or_default(json, "filter.elementwise_filter_type", &
         filter_type, "nonBoyd")
    this%filter_type = trim(filter_type)

    call json_get_or_default(json, "filter.direction", &
         direction, "rst")
    this%direction = trim(direction)
    if (this%direction .ne. "rst" .and. &
        this%direction .ne. "rs" .and. &
        this%direction .ne. "rt" .and. &
        this%direction .ne. "st" .and. &
        this%direction .ne. "r" .and. &
        this%direction .ne. "s" .and. &
        this%direction .ne. "t") then
       call neko_type_error("The direction of the elementwise " // &
            "filter", this%direction, KNOWN_DIRECTIONS)
    end if


    if (json%valid_path('filter.transfer_function')) then
       call json_get(json, 'filter.transfer_function', transfer)
       if (size(transfer) .eq. coef%dof%xh%lx) then
          this%transfer = transfer
       else
          call neko_error("The transfer function of the elementwise " // &
               "filter must correspond the order of the polynomial")
       end if
       call this%build_1d()
    end if

  end subroutine elementwise_filter_init_from_json

  !> Actual Constructor.
  !! @param nx number of points in an elements in one direction.
  subroutine elementwise_filter_init_from_components(this, nx, coef)
    class(elementwise_filter_t), intent(inout) :: this
    integer, intent(in) :: nx
    type(coef_t), intent(in), target, optional :: coef
    integer :: i

    this%nx = nx
    this%nt = nx ! initialize as if nothing is filtered yet
    
    ! Force the pointing of coef 
    ! since init_from_components is sometimes called alone
    if (present(coef)) then
       this%coef => coef
    end if

    allocate(this%fh(nx, nx))
    allocate(this%fht(nx, nx))
    allocate(this%ident(nx, nx))
    allocate(this%transfer(nx))

    call rzero(this%fh, nx*nx)
    call rzero(this%fht, nx*nx)
    call rzero(this%ident, nx*nx)
    call rone(this%transfer, nx) ! initialize as if nothing is filtered yet

    do i = 1, nx
       this%ident(i, i) = 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%fh, this%fh_d, this%nx * this%nx)
       call device_map(this%fht, this%fht_d, this%nx * this%nx)
       call device_map(this%ident, this%ident_d, this%nx * this%nx)
       call device_cfill(this%fh_d, 0.0_rp, this%nx * this%nx)
       call device_cfill(this%fht_d, 0.0_rp, this%nx * this%nx)
       call device_memcpy(this%ident, this%ident_d, &
            this%nx * this%nx, HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine elementwise_filter_init_from_components

  !> Destructor.
  subroutine elementwise_filter_free(this)
    class(elementwise_filter_t), intent(inout) :: this

    if (allocated(this%filter_type)) then
       deallocate(this%filter_type)
    end if

    if (allocated(this%fh)) then
       deallocate(this%fh)
    end if

    if (allocated(this%fht)) then
       deallocate(this%fht)
    end if

    if (allocated(this%transfer)) then
       deallocate(this%transfer)
    end if

    if (c_associated(this%fh_d)) then
       call device_free(this%fh_d)
    end if

    if (c_associated(this%fht_d)) then
       call device_free(this%fht_d)
    end if

    if (c_associated(this%ident_d)) then
       call device_free(this%ident_d)
    end if

    this%filter_type = ""
    this%nx = 0
    this%nt = 0

    call this%free_base()

  end subroutine elementwise_filter_free

  !> Build the 1d filter for an element.
  subroutine build_1d(this)
    class(elementwise_filter_t), intent(inout) :: this

    call build_1d_cpu(this%fh, this%fht, this%transfer, &
         this%nx, this%filter_type)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%fh, this%fh_d, &
            this%nx * this%nx, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%fht, this%fht_d, &
            this%nx * this%nx, HOST_TO_DEVICE, sync = .false.)
    end if

  end subroutine build_1d

  !> Filter a 3D field.
  subroutine elementwise_field_filter_3d(this, F_out, F_in)
    class(elementwise_filter_t), intent(inout) :: this
    type(field_t), intent(inout) :: F_out
    type(field_t), intent(in) :: F_in

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(F_out%x_d, F_out%dof%size())
    else
       call rzero(F_out%x, F_out%dof%size())
    end if

    ! F_out = fh x fh x fh x F_in
    if (this%direction .eq. "rst") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%fh, this%fht, &
            this%fht, this%coef%msh%nelv)
    else if (this%direction .eq. "rs") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%fh, this%fht, &
            this%ident, this%coef%msh%nelv)
    else if (this%direction .eq. "rt") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%fh, this%ident, &
            this%fht, this%coef%msh%nelv)
    else if (this%direction .eq. "st") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%ident, this%fht, &
            this%fht, this%coef%msh%nelv)
    else if (this%direction .eq. "r") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%fh, this%ident, &
            this%ident, this%coef%msh%nelv)
    else if (this%direction .eq. "s") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%ident, this%fht, &
            this%ident, this%coef%msh%nelv)
    else if (this%direction .eq. "t") then
       call tnsr3d(F_out%x, this%nx, F_in%x, this%nx, this%ident, this%ident, &
            this%fht, this%coef%msh%nelv)
    end if

  end subroutine elementwise_field_filter_3d

  !> Build the 1d filter for an element on the CPU.
  !> Suppose field x is filtered into x_hat by x_hat = fh*x.
  !! @param fh The 1D filter operator.
  !! @param fht The transpose of fh.
  !! @param trnfr The transfer function containing weights for different modes.
  !! @param nx number of points, dimension of x.
  !! @param filter_type
  subroutine build_1d_cpu(fh, fht, transfer, nx, filter_type)
    integer, intent(in) :: nx
    real(kind=rp), intent(inout) :: fh(nx, nx), fht(nx, nx)
    real(kind=rp), intent(in) :: transfer(nx)
    real(kind=rp) :: diag(nx, nx), rmult(nx), Lj(nx), zpts(nx)
    type(matrix_t) :: phi, pht
    integer :: n, i, j, k
    real(kind=rp) :: z
    character(len=*), intent(in) :: filter_type

    call phi%init(nx, nx)
    call pht%init(nx, nx)

    call zwgll(zpts, rmult, nx)

    n = nx-1
    do j = 1, nx
       z = zpts(j)
       call legendre_poly(Lj, z, n)
       select case (filter_type)
       case("Boyd")
          pht%x(1,j) = Lj(1)
          pht%x(2,j) = Lj(2)
          do k=3,nx
             pht%x(k,j) = Lj(k)-Lj(k-2)
          end do
       case("nonBoyd")
          pht%x(:,j) = Lj
       end select
    end do

    call trsp(phi%x, nx, pht%x, nx)
    pht%x = phi%x

    call pht%inverse(0) ! "0" for cpu implementation

    diag = 0.0_rp

    do i=1,nx
       diag(i,i) = transfer(i)
    end do

    call mxm (diag, nx, pht%x, nx, fh, nx) !          -1
    call mxm (phi%x, nx, fh, nx, pht%x, nx) !     V D V

    call copy (fh, pht%x, nx*nx)
    call trsp (fht, nx, fh, nx)

  end subroutine build_1d_cpu

end module elementwise_filter
