!> Defines a sampler 
module sampler
  use num_types
  use output
  use comm
  implicit none

  !> Pointer to an arbitrary output
  type, private :: outp_t
     class(output_t), pointer :: outp
  end type outp_t

  !> Sampler
  type :: sampler_t
     type(outp_t), allocatable :: output_list(:)
     integer :: n              
     integer :: size
     integer :: nsample
     real(kind=dp) :: freq
     real(kind=dp) :: T
   contains
     procedure, pass(this) :: init => sampler_init
     procedure, pass(this) :: free => sampler_free
     procedure, pass(this) :: add => sampler_add
     procedure, pass(this) :: sample => sampler_sample
  end type sampler_t

contains

  !> Initialize a sampler
  subroutine sampler_init(this, nsamp, T_end, size)
    class(sampler_t), intent(inout) :: this
    integer, intent(inout) :: nsamp
    real(kind=dp) :: T_end
    integer, intent(inout), optional :: size
    integer :: n, i

    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%output_list(n))

    do i = 1, n
       this%output_list(i)%outp => null()
    end do

    this%n = 0
    this%size = n

    this%nsample = 0
    this%freq = ( dble(nsamp) / T_end )
    this%T = 1d0 / this%freq
    
  end subroutine sampler_init

  !> Deallocate a sampler
  subroutine sampler_free(this)
    class(sampler_t), intent(inout) :: this

    if (allocated(this%output_list)) then
       deallocate(this%output_list)       
    end if

    this%n = 0
    this%size = 0

  end subroutine sampler_free

  !> Add an output @a out to the sampler
  subroutine sampler_add(this, out)
    class(sampler_t), intent(inout) :: this
    class(output_t), intent(inout), target :: out
    type(outp_t), allocatable :: tmp(:)

    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%output_list
       call move_alloc(tmp, this%output_list)
    end if

    this%n = this%n + 1
    this%output_list(this%n)%outp => out
    
  end subroutine sampler_add

  !> Sample all outputs in the sampler
  subroutine sampler_sample(this, t)
    class(sampler_t), intent(inout) :: this
    real(kind=dp), intent(in) :: t
    integer :: i

    if (t .ge. (this%nsample * this%T)) then

       if (pe_rank .eq. 0) then
          write(*,'(1x,a23,1x,e15.7)') 'Sampling fields at time:', t
       end if
       
       do i = 1, this%n
          call this%output_list(i)%outp%sample(t)
       end do
       this%nsample = this%nsample + 1
    end if
    
  end subroutine sampler_sample
  
end module sampler