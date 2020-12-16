!> Simulation driver
module simulation
  use case
  use abbdf
  implicit none
  private

  public :: neko_solve
  
contains

  !> Main driver to solve a case @a C
  subroutine neko_solve(C)
    type(case_t), intent(inout) :: C
    real(kind=dp) :: t
    integer :: i, tstep
    type(file_t) :: fileout

    t = 0d0
    fileout = file_t("oufluid.fld")
    do i = 1, C%params%nsteps
       tstep = i
       call simulation_settime(t, C%params%dt, C%ab_bdf, C%tlag, C%dtlag, tstep)
       call C%fluid%step(t, tstep, C%ab_bdf)
       if (pe_rank .eq. 0) write(*,*) 'Save'
       call fileout%write(C%fluid)
       !> @todo Add call to sampler
    end do
    
  end subroutine neko_solve

  subroutine simulation_settime(t, dt, ab_bdf, tlag, dtlag, step)
    real(kind=dp), intent(inout) :: t
    real(kind=dp), intent(in) :: dt
    type(abbdf_t), intent(inout) :: ab_bdf
    real(kind=dp), dimension(10) :: tlag
    real(kind=dp), dimension(10) :: dtlag
    integer, intent(in) :: step
    integer :: i

    do i = 10, 2, -1
       tlag(i) = tlag(i-1)
       dtlag(i) = dtlag(i-1)
    end do

    dtlag(1) = dt
    if (step .eq. 1) then
       dtlag(2) = dt
       tlag(2) = t
    end if

    t = t + dt

    call ab_bdf%set_bd(dtlag)
    call ab_bdf%set_abbd(dtlag)
    
  end subroutine simulation_settime

end module simulation

