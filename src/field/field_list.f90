module field_list
  use field, only : field_ptr_t, field_t
  use iso_c_binding, only : c_ptr
  implicit none
  private

  !> field_list_t, To be able to group fields together
  type, public :: field_list_t
     type(field_ptr_t), allocatable :: items(:)
   contains
     !> Constructor. Allocates array and pointers.
     procedure, pass(this) :: init => field_list_init
     !> Destructor.
     procedure, pass(this) :: free => field_list_free
     !> Append a field to the list.
     procedure, pass(this) :: append => field_list_append
     !> Get an item pointer by array index
     procedure, pass(this) :: get => field_list_get
     !> Point item at given index.
     generic :: set => set_to_ptr, set_to_field_ptr
     procedure, pass(this) :: set_to_ptr => field_list_set_to_ptr
     procedure, pass(this) :: set_to_field_ptr => field_list_set_to_field_ptr
     procedure, pass(this) :: set_to_field => field_list_set_to_field

     !> Get device pointer for a given index
     procedure, pass(this) :: x_d => field_list_x_d

     !> Get number of items in the list.
     procedure, pass(this) :: size => field_list_size
  end type field_list_t

contains
  !> Constructor. Just allocates the array.
  !! @param size The size of the list to preallocate
  subroutine field_list_init(this, size)
    class(field_list_t), intent(inout) :: this
    integer, intent(in) :: size
    integer :: i

    call this%free()

    allocate(this%items(size))
  end subroutine field_list_init

  !> Get number of items in the list.
  pure function field_list_size(this) result(n)
    class(field_list_t), intent(in) :: this
    integer :: n
    n = size(this%items)
  end function field_list_size

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function field_list_get(this, i) result(f)
    class(field_list_t), intent(inout) :: this
    type(field_t), pointer :: f
    integer :: i
    f = this%items(i)%ptr
  end function field_list_get

  !> Append a field to the list.
  !! @param f The field to append.
  subroutine field_list_append(this, f)
    class(field_list_t), intent(inout) :: this
    class(field_t), intent(in), target :: f
    type(field_ptr_t), allocatable :: tmp(:)
    integer :: len

    len = size(this%items)

    allocate(tmp(len+1))
    tmp(1:len) = this%items
    call move_alloc(tmp, this%items)
    this%items(len+1)%ptr => f

  end subroutine field_list_append

  !> Destructor.
  subroutine field_list_free(this)
    class(field_list_t), intent(inout) :: this
    integer :: i, n_fields

    if (allocated(this%items)) then
       n_fields = this%size()
       do i=1, n_fields
          call this%items(i)%ptr%free()
          nullify(this%items(i)%ptr)
       end do
       deallocate(this%items)
    end if

  end subroutine field_list_free

  !> Get device pointer for a given index
  !! @param i The index of the item.
  function field_list_x_d(this, i) result(ptr)
    class(field_list_t), intent(inout) :: this
    type(c_ptr), pointer :: ptr
    integer :: i
    ptr = this%items(i)%ptr%x_d
  end function field_list_x_d

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr A field pointer to point the item to.
  subroutine field_list_set_to_ptr(this, i, ptr)
    class(field_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_t), pointer, intent(in) :: ptr

    this%items(i)%ptr => ptr
  end subroutine field_list_set_to_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr An encapsulated field pointer to point the item to.
  subroutine field_list_set_to_field_ptr(this, i, ptr)
    class(field_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_ptr_t), intent(in) :: ptr

    this%items(i)%ptr => ptr%ptr
  end subroutine field_list_set_to_field_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param field A field to point the item to.
  subroutine field_list_set_to_field(this, i, fld)
    class(field_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_t), target, intent(in) :: fld

    this%items(i)%ptr => fld
  end subroutine field_list_set_to_field



end module field_list
