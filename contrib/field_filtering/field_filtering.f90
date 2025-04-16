!> Program to filter a field
program field_filtering
  use neko
  use filter, only : filter_t
  use Najafi_Yazdi_filter, only : Najafi_Yazdi_filter_t
  use PDE_filter, only : PDE_filter_t
  use AHO_procedure, only : AHO_procedure_t
  use scratch_registry, only : scratch_registry_t, neko_scratch_registry
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field_fname, output_fname, mesh_fname
  type(file_t) :: field_file, output_file, mesh_file
  type(fld_file_data_t) :: field_data
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(coef_t), target :: coef
  type(dofmap_t) :: dof
  type(gs_t) :: gs_h
  type(vector_ptr_t), allocatable :: fields(:)
  type(field_t) :: field_in, field_out
  integer :: argc, i, lx, j, file_precision
  logical :: dp_precision
!   type(PDE_filter_t) :: filter
  ! type(Najafi_Yazdi_filter_t) :: filter
  type(AHO_procedure_t) :: filter
  real(kind=rp) :: t_start, t_end, t_elapsed
  t_elapsed = 0.0_rp

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

  neko_scratch_registry = scratch_registry_t(dof, 10, 10)

  !! Initialize the PDE filter
  ! for order 9
! for PDE filter
!   filter%r = 0.031746031746031744 !! Transfer function 0.5 at average GLL
!   filter%r = 0.047222559333253436 !! Transfer function 0.5 at max GLL
  
  ! filter%alpha = 0.027267807205288697 !! Transfer function 0.5 at average GLL
! !   filter%alpha = 0.04433613524144137 !! Transfer function 0.5 at max GLL
!   filter%alpha = 0.016256605483164713 !! Transfer function 0.5 at double min GLL
! !   filter%beta = 0.01149515597622018 !! zero at the minimal GLL spacing
  filter%beta = 0.01149515597622018 !! zero at the minimal GLL spacing
! !   filter%beta = 0.0
  filter%AD_order = 3
  filter%gamma = 2.212
  filter%damp_order = 8
  allocate(Najafi_Yazdi_filter_t::filter%base_filter)
  
  select type(f => filter%base_filter)
  type is (Najafi_Yazdi_filter_t)
    f%alpha = 0.027267807205288697 !! Transfer function 0.5 at average GLL
    ! f%alpha = 0.04433613524144137 !! Transfer function 0.5 at max GLL
   !  f%alpha = 0.016256605483164713 !! Transfer function 0.5 at double min GLL
    !   filter%base_filter%beta = 0.01149515597622018 !! zero at the minimal GLL spacing
    f%beta = 0.01149515597622018 !! zero at the minimal GLL spacing
    !   filter%base_filter%beta = 0.0
  end select

  ! general info for filter
  ! filter%abstol_filt = 1e-10
  ! filter%ksp_max_iter = 800
  ! filter%ksp_solver = "cg"
  ! filter%precon_type_filt = 'jacobi'
  select type(f => filter%base_filter)
  type is (Najafi_Yazdi_filter_t)
     f%abstol_filt = 1e-10
     f%ksp_max_iter = 800
     f%ksp_solver = "cg"
     f%precon_type_filt = 'jacobi'
     f%coef => coef
     call f%init_from_attributes(coef)
  end select
  filter%coef => coef
  call filter%init_from_attributes(coef)

  ! interpolate field for t>0
  allocate(fields(field_data%size()))
  do i = 1, field_data%meta_nsamples
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i
     if (i .gt. 1) then
        call field_file%read(field_data)
     end if
     call field_data%get_list(fields,field_data%size())
     if (pe_rank .eq. 0) call cpu_time(t_start)
     do j = 1, field_data%size()
        !! Apply filter to field fields(j)%ptr%x
        call copy(field_in%x, fields(j)%ptr%x, field_in%dof%size())
        call filter%apply(field_out, field_in)
        call copy(fields(j)%ptr%x, field_out%x, field_in%dof%size())
     end do
     if (pe_rank .eq. 0) then
        call cpu_time(t_end)
        t_elapsed = t_elapsed + t_end - t_start
        write(*,*) "elapsed time", t_end - t_start 
     end if
     ! output
     call output_file%write(field_data, field_data%time)
  end do
  if (pe_rank .eq. 0) write(*,*) "Total elapsed time for filtering: ", t_elapsed 


  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program field_filtering
