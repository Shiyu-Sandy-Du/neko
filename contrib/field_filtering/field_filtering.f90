!> Program to filter a field
program field_filtering
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, output_file
  type(fld_file_data_t) :: field_data
  type(space_t) :: Xh
  type(vector_ptr_t), allocatable :: fields(:)
  type(field_t), allocatable :: field_in, field_out
  integer :: argc, i, lx, j, file_precision
  logical :: dp_precision

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./field_filtering field.fld outfield.fld precision'
        write(*,*) 'Example command: ./field_filtering fieldblabla.fld outfield.fld .true.'
        write(*,*) 'Filter a series of field'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, fmt='(A)') field_fname
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') output_fname
  call get_command_argument(3, inputchar)
  read(inputchar, *) dp_precision

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  field_file = file_t(trim(field_fname),precision=file_precision)

  call field_data%init()
  if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
  call field_file%read(field_data)

  lx = field_data%lx

  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)

  ! interpolate field for t>0
  do i = 1, field_data%meta_nsamples
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i
     call field_file%read(field_data)
     call field_data%get_list(fields,field_data%size())
     do j = 1, field_data%size()
        !! Apply filter to field fields(j)%ptr%x
     end do
     ! output for t>0
     call output_file%write(field_data, field_data%time)
  end do


  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program field_filtering
