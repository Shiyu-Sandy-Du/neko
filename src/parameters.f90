!> NEKO parameters
module parameters
  use num_types
  implicit none  

  type param_t
     integer :: nsamples        !< Number of samples
     logical :: output_bdry     !< Output boundary markings
     real(kind=dp) :: dt        !< time-step size               
     real(kind=dp) :: T_end     !< Final time
     real(kind=dp) :: rho       !< Density \f$ \rho \f$
     real(kind=dp) :: mu        !< Dynamic viscosity \f$ \mu \f$
     real(kind=dp) :: Re        !< Reynolds number
     real(kind=dp), dimension(3) :: uinf !< Free-stream velocity \f$ u_\infty \f$
  end type param_t

  type param_io_t
     type(param_t) p
   contains
     procedure  :: param_read
     generic :: read(formatted) => param_read
  end type param_io_t

  interface write(formatted)
     module procedure :: param_write
  end interface write(formatted)
  
contains

  subroutine param_read(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(inout) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    integer :: nsamples = 0
    logical :: output_bdry = .false.
    real(kind=dp) :: dt = 0d0
    real(kind=dp) :: T_end = 0d0
    real(kind=dp) :: rho = 1d0
    real(kind=dp) :: mu = 1d0
    real(kind=dp) :: Re = 1d0
    real(kind=dp), dimension(3) :: uinf = (/ 0d0, 0d0, 0d0 /)

    namelist /NEKO_PARAMETERS/ nsamples, output_bdry, dt, T_end, rho, mu, &
         Re, uinf

    read(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)

    param%p%output_bdry = output_bdry
    param%p%nsamples = nsamples
    param%p%dt = dt
    param%p%T_end = T_end
    param%p%rho = rho 
    param%p%mu = mu
    param%p%Re = Re
    param%p%uinf = uinf


  end subroutine param_read

  subroutine param_write(param, unit, iotype, v_list, iostat, iomsg)
    class(param_io_t), intent(in) ::  param
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer(kind=4), intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    real(kind=dp) :: dt, T_End, rho, mu, Re
    real(kind=dp), dimension(3) :: uinf
    logical :: output_bdry
    integer :: nsamples
    namelist /NEKO_PARAMETERS/ nsamples, output_bdry, dt, T_end, rho, mu, &
         Re, uinf

    nsamples = param%p%nsamples
    output_bdry = param%p%output_bdry
    dt = param%p%dt
    T_end = param%p%T_end
    rho = param%p%rho  
    mu = param%p%mu
    Re = param%p%Re
    uinf = param%p%uinf

    
    write(unit, nml=NEKO_PARAMETERS, iostat=iostat, iomsg=iomsg)

        
  end subroutine param_write

  
end module parameters

