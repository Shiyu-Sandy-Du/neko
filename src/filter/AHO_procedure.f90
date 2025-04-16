! Copyright (c) 2025, The Neko Authors
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
!> An approximate high order (AHO) procedure for 2nd order differential filters
module AHO_procedure
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get_or_default, json_get, json_extract_object
  use field_registry, only: neko_field_registry
  use field, only: field_t
  use coefs, only: coef_t
  use ax_product, only: ax_t, ax_helm_factory
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory, &
       krylov_solver_destroy
  use precon, only: pc_t, precon_factory, precon_destroy
  use bc_list, only : bc_list_t
  use neumann, only: neumann_t
  use profiler, only: profiler_start_region, profiler_end_region
  use gather_scatter, only: gs_t, GS_OP_ADD
  use pnpn_residual, only: pnpn_prs_res_t
  use mesh, only: mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use field_registry, only: neko_field_registry
  use filter, only: filter_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_copy, field_add2, field_sub2, field_add2s2
  use coefs, only: coef_t
  use logger, only: neko_log, LOG_SIZE
  use neko_config, only: NEKO_BCKND_DEVICE
  use dofmap, only: dofmap_t
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use sx_jacobi, only: sx_jacobi_t
  use hsmg, only: hsmg_t
  use utils, only : concat_string_array, neko_error
  use math, only : invcol2, col2
  use device_math, only: device_cfill, device_subcol3, &
                         device_cmult, device_invcol2, &
                         device_col2
  use elementwise_filter, only : elementwise_filter_t
  use PDE_filter, only : PDE_filter_t
  use Najafi_Yazdi_filter, only : Najafi_Yazdi_filter_t
  implicit none
  private

  !> An approximate high order (AHO) procedure based on pde filter
  type, public, extends(filter_t) :: AHO_procedure_t
     !> Base filter
     class(filter_t), allocatable :: base_filter 
     !> Ax
     class(ax_t), allocatable :: Ax_damp
     ! Inputs from the user
     !> Order of the approximate deconvolution (AD)-like procedure
     integer :: AD_order
     !> Magnifying coefficient for the higher-order term in AD
     real(kind=rp) :: gamma
     !> Order of the damp procedure 
     integer :: damp_order
     !> The parameter determining the vanishing wavenumber
     real(kind=rp) :: beta


   contains
     !> Constructor from json.
     procedure, pass(this) :: init => AHO_procedure_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          AHO_procedure_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => AHO_procedure_free
     !> Apply the filter
     procedure, pass(this) :: apply => AHO_procedure_apply
  end type AHO_procedure_t

contains

  !> Constructor from json.
  subroutine AHO_procedure_init_from_json(this, json, coef)
    class(AHO_procedure_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(json_file) :: json_base_filter
    type(coef_t), intent(in) :: coef
    character(len=:), allocatable :: base_filter_type
    character(len=:), allocatable :: type_string
    character(len=20) :: FILTER_KNOWN_TYPES(4) = [character(len=20) :: &
     "elementwise", &
     "PDE", &
     "Yazdi", &
     "AHO"]

    ! user parameters
    call json_get(json, "filter.AD_order", this%AD_order)
    call json_get(json, "filter.gamma", this%gamma)
    call json_get(json, "filter.damp_order", this%damp_order)
    call json_get(json, "filter.beta", this%beta)

    call this%init_base(json, coef)

    ! initialize the base filter
    call json_extract_object(json, "filter.base_filter", json_base_filter)
    call json_get(json_base_filter, 'filter.type', base_filter_type)

    if (allocated(this%base_filter)) then
       deallocate(this%base_filter)
    else if (trim(base_filter_type) .eq. 'elementwise') then
       allocate(elementwise_filter_t::this%base_filter)
    else if (trim(base_filter_type) .eq. 'PDE') then
       allocate(pde_filter_t::this%base_filter)
    else if (trim(base_filter_type) .eq. 'Najafi_Yazdi') then
       allocate(Najafi_Yazdi_filter_t::this%base_filter)
    else if (trim(base_filter_type) .eq. 'AHO') then
       allocate(AHO_procedure_t::this%base_filter)
    else
       type_string =  concat_string_array(FILTER_KNOWN_TYPES, &
            NEW_LINE('A') // "-  ", .true.)
       call neko_error("Unknown filter type: " &
                       // trim(base_filter_type) // ".  Known types are: " &
                       // type_string)
       stop

    end if
    call this%base_filter%init(json_base_filter, coef)
    ! initialize the rest used in the filter
    call AHO_procedure_init_from_attributes(this, coef)

  end subroutine AHO_procedure_init_from_json

  !> Actual constructor.
  subroutine AHO_procedure_init_from_attributes(this, coef)
    class(AHO_procedure_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer :: n

    n = this%coef%dof%size()

    ! Setup backend dependent Ax routines
    call ax_helm_factory(this%Ax_damp, full_formulation = .false.)

  end subroutine AHO_procedure_init_from_attributes

  !> Destructor.
  subroutine AHO_procedure_free(this)
    class(AHO_procedure_t), intent(inout) :: this
    
    if (allocated(this%base_filter)) then
       call this%base_filter%free
       deallocate(this%base_filter)
    end if

    if (allocated(this%Ax_damp)) then
       deallocate(this%Ax_damp)
    end if

    call this%free_base()

  end subroutine AHO_procedure_free

  !> Apply the filter
  !! @param F_out filtered field
  !! @param F_in unfiltered field
  subroutine AHO_procedure_apply(this, F_out, F_in)
    class(AHO_procedure_t), intent(inout) :: this
    type(field_t), intent(in) :: F_in
    type(field_t), intent(inout) :: F_out
    integer :: n, i
    type(field_t) :: tmp_field, tmp_field2

    n = this%coef%dof%size()
    call tmp_field%init(this%coef%dof)
    call tmp_field2%init(this%coef%dof)
    call field_copy(tmp_field, F_in)

    !! Step 1: Pre-apply the damp procedure to smoothen the field
    ! set up Helmholtz operators for \rho + beta^2 \nabla^2 \rho
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%coef%h1_d, -this%beta**2, n)
       call device_cfill(this%coef%h2_d, 1.0_rp, n)
    else
       do i = 1, n
          ! h1 is already negative in its definition
          this%coef%h1(i,1,1,1) = -this%beta**2
          ! ax_helm includes the mass matrix in h2
          this%coef%h2(i,1,1,1) = 1.0_rp
       end do
    end if
    this%coef%ifh2 = .true.

    do i = 1, this%damp_order
       call this%Ax_damp%compute(F_out%x, tmp_field%x, &
       this%coef, this%coef%msh, this%coef%Xh)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(F_out%x_d, this%coef%B_d, n)
       else
          call invcol2(F_out%x, this%coef%B, n)
       end if
       ! gather scatter
       call this%coef%gs_h%op(F_out, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(F_out%x_d, this%coef%mult_d, n)
       else
          call col2(F_out%x, this%coef%mult, n)
       end if
       call field_copy(tmp_field, F_out)
    end do    

    !! Step 2: Apply the base filter
    call this%base_filter%apply(F_out, tmp_field)

    !! Step 3: Apply the approximate devoncolution procedure
    call field_copy(tmp_field, F_out)
    do i = 1, this%AD_order
       call this%base_filter%apply(tmp_field2, tmp_field)
       call field_sub2(tmp_field, tmp_field2)
       call field_add2s2(F_out, tmp_field, this%gamma)
    end do 
  end subroutine AHO_procedure_apply

end module AHO_procedure
