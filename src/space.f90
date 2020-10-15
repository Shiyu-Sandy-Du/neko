!> Defines a function space
module space
  use num_types
  use speclib
  use utils
  implicit none

  integer, parameter :: GL = 0, GLL = 1, GJ = 2

  type space_t
     integer :: lx              !< Polynomial dimension in x-direction
     integer :: ly              !< Polynomial dimension in y-direction
     integer :: lz              !< Polynomial dimension in z-direction
     
     real(kind=dp), allocatable :: zg(:,:) !< Quadrature points

     real(kind=dp), allocatable :: wx(:)   !< Quadrature weights
     real(kind=dp), allocatable :: wy(:)   !< Quadrature weights
     real(kind=dp), allocatable :: wz(:)   !< Quadrature weights

     real(kind=dp), allocatable :: w3(:,:,:)

     real(kind=dp), allocatable :: dx(:,:) !< Derivative operator
     real(kind=dp), allocatable :: dy(:,:) !< Derivative operator
     real(kind=dp), allocatable :: dz(:,:) !< Derivative operator

     real(kind=dp), allocatable :: dxt(:,:) !< Derivative operator
     real(kind=dp), allocatable :: dyt(:,:) !< Derivative operator
     real(kind=dp), allocatable :: dzt(:,:) !< Derivative operator

     !> @todo Store gll points etc in the space
  end type space_t

  interface operator(.eq.)
     module procedure space_eq
  end interface operator(.eq.)

  interface operator(.ne.)
     module procedure space_eq
  end interface operator(.ne.)
  
contains

  !> Initialize a function space @a s with given polynomial dimensions
  subroutine space_init(s, t, lx, ly, lz)
    type(space_t), intent(inout) :: s
    integer, intent(in) :: t    !< Quadrature type
    integer, intent(in) :: lx   !< Polynomial dimension in x-direction
    integer, intent(in) :: ly   !< Polynomial dimension in y-direction
    integer, intent(in) :: lz   !< Polynomial dimension in z-direction
    integer ::ix, iy, iz

    call space_free(s)

    s%lx = lx
    s%ly = ly
    s%lz = lz

    allocate(s%zg(lx, 3))

    allocate(s%wx(lx))
    allocate(s%wy(ly))
    allocate(s%wz(lz))

    allocate(s%w3(lx, ly, lz))

    allocate(s%dx(lx, lx))
    allocate(s%dy(ly, ly))
    allocate(s%dz(lz, lz))

    allocate(s%dxt(lx, lx))
    allocate(s%dyt(ly, ly))
    allocate(s%dzt(lz, lz))
    
    !>@todo add 2d case
    if (t .eq. GLL) then
       call zwgll(s%zg(1,1), s%wx, lx)
       call zwgll(s%zg(1,2), s%wy, ly)
       call zwgll(s%zg(1,3), s%wz, lz)
    else if (t .eq. GL) then
       call zwgl(s%zg(1,1), s%wx, lx)
       call zwgl(s%zg(1,2), s%wy, ly)
       call zwgl(s%zg(1,3), s%wz, lz)
    else
       call neko_error("Invalid quadrature rule")
    end if

    do iz = 1, lz
       do iy = 1, ly
          do ix = 1, lx
             s%w3(ix, iy, iz) = s%wx(ix) * s%wy(iy) * s%wz(iz)
          end do
       end do
    end do

    call dgll(s%dx, s%dxt, s%zg(1,1), lx, lx)
    call dgll(s%dy, s%dyt, s%zg(1,2), ly, ly)
    call dgll(s%dz, s%dzt, s%zg(1,3), lz, lz)
  end subroutine space_init
   
  !> Deallocate a space @a s
  subroutine space_free(s)
    type(space_t), intent(inout) :: s
    
    if (allocated(s%zg)) then
       deallocate(s%zg)
    end if

    if (allocated(s%wx)) then
       deallocate(s%wx)
    end if

    if (allocated(s%wy)) then
       deallocate(s%wy)
    end if

    if (allocated(s%wz)) then
       deallocate(s%wz)
    end if

    if (allocated(s%w3)) then
       deallocate(s%w3)
    end if

    if (allocated(s%dx)) then
       deallocate(s%dx)
    end if

    if (allocated(s%dy)) then
       deallocate(s%dy)
    end if

    if (allocated(s%dz)) then
       deallocate(s%dz)
    end if

    if (allocated(s%dxt)) then
       deallocate(s%dxt)
    end if

    if (allocated(s%dyt)) then
       deallocate(s%dyt)
    end if

    if (allocated(s%dzt)) then
       deallocate(s%dzt)
    end if

  end subroutine space_free

  !> Check if \f$ X_h = Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_eq(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Xh%lx) .and. &
         (Xh%ly .eq. Xh%ly) .and. &
         (Xh%lz .eq. Xh%lz) ) then
       res = .true.
    else
       res = .false.
    end if
    
  end function space_eq

  !> Check if \f$ X_h \ne Y_H \f$
  !! @note this only checks the polynomial dimensions
  pure function space_ne(Xh, Yh) result(res)
    type(space_t), intent(in) :: Xh
    type(space_t), intent(in) :: Yh
    logical :: res

    if ( (Xh%lx .eq. Xh%lx) .and. &
         (Xh%ly .eq. Xh%ly) .and. &
         (Xh%lz .eq. Xh%lz) ) then
       res = .false.
    else
       res = .true.
    end if
    
  end function space_ne


end module space
