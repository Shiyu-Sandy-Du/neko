!> Program to filter a field
program field_diff
  use neko
  use PDE_filter, only : PDE_filter_t
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field1_fname, field2_fname, &
                                   output_fname
  type(file_t) :: field1_file, field2_file, output_file
  type(fld_file_data_t) :: field1_data, field2_data
  type(vector_ptr_t), allocatable :: fields1(:), fields2(:)
  integer :: argc, i, lx, j, file_precision
  logical :: dp_precision

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./field_diff field1.fld field2.fld &
                     &outfield.fld precision'
        write(*,*) 'Example command: ./field_diff field1.fld &
                     &field2.fld outfield.fld .true.'
        write(*,*) 'Make the difference of two series of field'
     end if
     stop
  end if

  call neko_init
  
  call get_command_argument(1, inputchar)
  read(inputchar, fmt='(A)') field1_fname
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') field2_fname
  call get_command_argument(3, inputchar)
  read(inputchar, fmt='(A)') output_fname
  call get_command_argument(4, inputchar)
  read(inputchar, *) dp_precision

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  field1_file = file_t(trim(field1_fname),precision=file_precision)
  field2_file = file_t(trim(field2_fname),precision=file_precision)
  output_file = file_t(trim(output_fname),precision=file_precision)

  call field1_data%init()
  call field2_data%init()

  call field1_file%read(field1_data)
  call field2_file%read(field2_data)

  ! interpolate field for t>0
  allocate(fields1(field1_data%size()))
  allocate(fields2(field1_data%size()))
  do i = 1, field1_data%meta_nsamples
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i
     if (i .gt. 1) then
        call field1_file%read(field1_data)
        call field2_file%read(field2_data)
     end if
     call field1_data%get_list(fields1,field1_data%size())
     call field2_data%get_list(fields2,field1_data%size())
     do j = 1, field1_data%size()
        call sub2(fields1(j)%ptr%x, fields2(j)%ptr%x, size(fields1(j)%ptr%x))
     end do
     ! output
     call output_file%write(field1_data, field1_data%time)
  end do


  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program field_diff
