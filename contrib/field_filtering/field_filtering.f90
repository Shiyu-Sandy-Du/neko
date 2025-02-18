!> Program to filter a field
program field_filtering
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field_fname, hom_dir, output_fname, mesh_fname
  type(file_t) :: field_file, output_file, mesh_file
  type(fld_file_data_t) :: field_data
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(coef_t) :: coef
  type(dofmap_t) :: dof
  type(gs_t) :: gs_h
  type(vector_ptr_t), allocatable :: fields(:)
  type(field_t) :: field_in, field_out
  integer :: argc, i, lx, j, file_precision
  logical :: dp_precision
  real(kind=rp) :: r

  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./field_filtering mesh.nmsh field.fld outfield.fld precision'
        write(*,*) 'Example command: ./field_filtering box.nmsh fieldblabla.fld outfield.fld .true.'
        write(*,*) 'Filter a series of field'
     end if
     stop
  end if

  call neko_init
  
  call get_command_argument(1, inputchar)
  read(inputchar, *) mesh_fname
  mesh_file = file_t(trim(mesh_fname))
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') field_fname
  call get_command_argument(3, inputchar)
  read(inputchar, fmt='(A)') output_fname
  call get_command_argument(4, inputchar)
  read(inputchar, *) dp_precision

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  field_file = file_t(trim(field_fname),precision=file_precision)
  output_file = file_t(trim(output_fname),precision=file_precision)

  call field_data%init()
  if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
  call field_file%read(field_data)

  !! Initiate the work field_t for the filtered and unfiltered field
  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)
  call mesh_file%read(msh)
  call dof%init(msh, Xh)
  call gs_h%init(dof)
  call coef%init(gs_h)
  call field_in%init(dof, "field_in")
  call field_out%init(dof, "field_out")

  ! interpolate field for t>0
  call field_file%read(field_data)
  allocate(fields(field_data%size()))
  do i = 1, field_data%meta_nsamples - 1
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i
     if (i .gt. 1) then
        call field_file%read(field_data)
     end if
     call field_data%get_list(fields,field_data%size())
     do j = 1, field_data%size()
        !! Apply filter to field fields(j)%ptr%x
        call copy(field_in%x, fields(j)%ptr%x, field_in%dof%size())
        call copy(field_out%x, field_in%x, field_in%dof%size())
        call copy(fields(j)%ptr%x, field_out%x, field_in%dof%size())
     end do
     ! output
     call output_file%write(field_data, field_data%time)
  end do


  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program field_filtering
