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
!> A simulation component that computes residual_viscosity
!! The values are stored in the field registry under the name 'residual_viscosity'

module residual_viscosity
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use scratch_registry, only : neko_scratch_registry
  use json_utils, only : json_get
  use field, only : field_t, field_ptr_t
  use rhs_maker, only : rhs_maker_bdf_t, rhs_maker_ext_t
  use advection, only : advection_t
  use time_scheme_controller, only : time_scheme_controller_t
  use coefs, only : coef_t
  use field_series, only : field_series_ptr_t
  use time_state, only : time_state_t
  use case, only : case_t
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use fluid_pnpn, only : fluid_pnpn_t
  use scalars, only : scalars_t
  use utils, only : neko_error
  use field_math, only : field_col3, field_copy, field_absval, field_rzero, &
                         field_cmult, field_sub2, field_col2, field_cadd2, &
                         field_invcol2
  use math, only : invcol2, col2, glsum, glsc2, glmax
  use device_math, only : device_invcol2
  use gather_scatter, only : GS_OP_ADD
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: residual_viscosity_t
     !> coefficient
     real(kind=rp) :: c_E
     !> X velocity component.
     type(field_t), pointer :: u
     type(field_t) :: u_var
     !> Y velocity component.
     type(field_t), pointer :: v
     type(field_t) :: v_var
     !> Z velocity component.
     type(field_t), pointer :: w
     type(field_t) :: w_var
     !> Velocity magnitude.
     type(field_t) :: vel_mag
     !> Scalar field
     integer :: n_scalars = 0
     type(field_ptr_t), allocatable :: s(:)
     type(field_t), allocatable :: s_var(:)
     !> coef
     type(coef_t), pointer :: coef
     !> Some field to be used
     type(field_t), allocatable :: abx1(:), abx2(:)
     type(field_t) :: h_np2
     real(kind=rp) :: volume_domain

     !> X residual_viscosity component.
     type(field_ptr_t), allocatable :: residual_viscosity(:)

     !> Residual.
     type(field_t), allocatable :: D(:)

     !> Output writer.
     type(field_writer_t) :: writer

     !> A pointer pointing to the time scheme
     type(fluid_pnpn_t), pointer :: fluid
     type(scalars_t), pointer :: scalars
     class(rhs_maker_bdf_t), pointer :: makebdf
     class(rhs_maker_ext_t), pointer :: makeext
     class(advection_t), pointer :: adv
     type(time_scheme_controller_t), pointer :: ext_bdf

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => residual_viscosity_init_from_json
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => residual_viscosity_init_common
     !> Destructor.
     procedure, pass(this) :: free => residual_viscosity_free
     !> Part of the residual viscosity computation before the time stepping
     procedure, pass(this) :: preprocess_ => residual_viscosity_preprocess
     !> Part of the residual viscosity computation after the time stepping
     procedure, pass(this) :: compute_ => residual_viscosity_compute
  end type residual_viscosity_t

contains

  !> Constructor from json.
  subroutine residual_viscosity_init_from_json(this, json, case)
    class(residual_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case

    call this%init_base(json, case)
    
    call json_get(json, "c_E", this%c_E)

    call this%init_common(json, case)
  end subroutine residual_viscosity_init_from_json

  !> Common part of constructors.
  subroutine residual_viscosity_init_common(this, json, case)
    class(residual_viscosity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20), allocatable :: fields(:)
    integer :: e, k
    real(kind=rp) :: volume_element

    if (allocated(case%scalars)) then
       this%scalars => case%scalars
       this%n_scalars = size(this%scalars%scalar_fields)
    else
       this%n_scalars = 0
    end if
    allocate(fields(3+this%n_scalars))
    fields(1) = 'res_visc_u'
    fields(2) = 'res_visc_v'
    fields(3) = 'res_visc_w'
    do k = 1, this%n_scalars
       write(fields(k+3), '(A,I0)') 'res_visc_s', k
    end do
    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to 	simulation_components/residual_viscosity.f90\the registry.
    call json%add("fields", fields)
    call this%writer%init(json, case)

    this%coef => case%fluid%c_Xh

    select type (f1 => case%fluid)
    type is (fluid_pnpn_t)
      this%fluid => f1
      this%makebdf => f1%makebdf
      this%makeext => f1%makeabf
      this%adv => f1%adv
      this%ext_bdf => f1%ext_bdf
    class default
      call neko_error("For fluid, residual &
      &viscosity currently only support pnpn scheme")
    end select

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    call this%u_var%init(u%dof)
    call this%v_var%init(v%dof)
    call this%w_var%init(w%dof)
    call this%vel_mag%init(u%dof)

    call this%h%init(this%u%dof)

    if (this%n_scalars .ne. 0) then
       allocate(this%s(this%n_scalars))
       allocate(this%s_var(this%n_scalars))
    end if

    allocate(this%D(3+this%n_scalars))
    allocate(this%abx1(3+this%n_scalars))
    allocate(this%abx2(3+this%n_scalars))
    
    allocate(this%residual_viscosity(3+this%n_scalars))

    do k = 1, 3+this%n_scalars
       this%residual_viscosity(k)%ptr => &
              neko_field_registry%get_field(fields(k))

       call this%D(k)%init(this%u%dof)
       call this%abx1(k)%init(this%u%dof)
       call this%abx2(k)%init(this%u%dof)
       
       if (k .le. this%n_scalars) then
          this%s(k)%ptr => this%scalars%scalar_fields(k)%s
          call this%s_var(k)%init(this%u%dof)
       end if
    end do

    do e = 1, this%coef%msh%nelv
       volume_element = 0.0_rp
       do k = 1, this%coef%Xh%lx * this%coef%Xh%ly * this%coef%Xh%lz
          volume_element = volume_element + this%coef%B(k, 1, 1, e)
       end do
       this%h_np2%x(:,:,:,e) = (volume_element**(1.0_rp/3.0_rp) &
                          / (this%coef%Xh%lx-1.0_rp)) &
                          **(this%ext_bdf%diffusion_time_order + 2)
    end do

    volume_domain = glsum(this%coef%B, u%dof%size())

  end subroutine residual_viscosity_init_common

  !> Destructor.
  subroutine residual_viscosity_free(this)
    class(residual_viscosity_t), intent(inout) :: this
    call this%free_base()
  end subroutine residual_viscosity_free

  !> Part of the residual_viscosity computation before the time stepping.
  !! @param time The time state.
  subroutine residual_viscosity_preprocess(this, time)
    class(residual_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, n

    !! estimate by the difference of the extrapolated and the solved advection
    associate(u => this%u, v => this%v, w => this%w, &
              D => this%D, &
              abx1 => this%abx1, abx2 => this%abx2, &
              coef => this%coef, &
              rho => this%fluid%rho, dt => time%dt, &
              adv => this%adv, &
              makeext => this%makeext, ext_bdf => this%ext_bdf, &
              Xh => this%coef%Xh)
  
    n = u%dof%size()
    call field_rzero(D(1))
    call field_rzero(D(2))
    call field_rzero(D(3))
    call adv%compute(u, v, w, D(1), D(2), D(3), Xh, coef, n)
    call makeext%compute_fluid(abx1(1), abx1(2), abx1(3), abx2(1), abx2(2), &
            abx2(3), D(1)%x, D(2)%x, D(3)%x, &
            rho%x(1,1,1,1), ext_bdf%advection_coeffs, n)

    end associate
    
    do i = 1, this%n_scalars
       !! estimate by the difference of the extrapolated and the solved advection
       associate(u => this%u, v => this%v, w => this%w, &
                 s => this%s(i)%ptr, D => this%D(i+3), &
                 abx1 => this%abx1(i+3), abx2 => this%abx2(i+3), &
                 coef => this%coef, &
                 rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
                 adv => this%adv, &
                 makeext => this%makeext, ext_bdf => this%ext_bdf, &
                 Xh => this%coef%Xh)
      
       n = s%dof%size()
       call field_rzero(D)
       call adv%compute_scalar(u, v, w, s, D, &
                Xh, coef, n)
       call makeext%compute_scalar(abx1, abx2, D%x, &
                rho%x(1,1,1,1), ext_bdf%advection_coeffs, n)

       end associate
       !  !! estimate by ds/dt + ui ds/dxi
      !  associate(s => this%s(i)%ptr, wa => this%wa(i), &
      !            coef => this%coef, &
      !            rho => this%scalars%scalar_fields(i)%rho, dt => time%dt, &
      !            makebdf => this%makebdf, &
      !            slag => this%slag(i)%ptr, ext_bdf => this%ext_bdf)
      
      !  n = s%dof%size()
      !  call field_rzero(wa)
      !  call makebdf%compute_scalar(slag, wa%x, s, coef%B, rho%x(1,1,1,1), &
      !          dt, ext_bdf%diffusion_coeffs, ext_bdf%ndiff, n)
      !  end associate
    end do

  end subroutine residual_viscosity_preprocess

  !> Part of the residual_viscosity computation after the time stepping.
  !! @param time The time state.
  subroutine residual_viscosity_compute(this, time)
    class(residual_viscosity_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_ptr_t) :: ta(3+this%n_scalars) ! temporal array
    real(kind=rp) :: u_avg, v_avg, w_avg
    real(kind=rp) :: s_avg(this%n_scalars)
    integer :: temp_indices(3+this%n_scalars)
    integer :: i, j, n
    real(kind=rp) :: scaling_uvw(3), scaling_s(this%n_scalars)

    do i = 1, 3+this%n_scalars
       call neko_scratch_registry%request_field(ta(i)%ptr, temp_indices(i))
    end do

    !! estimate by the difference of the extrapolated and the solved advection
    associate(ext_bdf => this%ext_bdf, &
              dt => time%dt, coef => this%coef, D => this%D, &
              adv => this%adv, gs => this%coef%gs_h, ta => ta, &
              u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
              residual_viscosity => this%residual_viscosity)

    n = u%dof%size()
    call field_rzero(ta(1)%ptr)
    call field_rzero(ta(2)%ptr)
    call field_rzero(ta(3)%ptr)
    call adv%compute(u, v, w, ta(1)%ptr, ta(2)%ptr, ta(3)%ptr, Xh, coef, n)
    call field_sub2(D(1), ta(1)%ptr, n)
    call field_sub2(D(2), ta(2)%ptr, n)
    call field_sub2(D(3), ta(3)%ptr, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_invcol2(D(1)%x_d, coef%B_d, n)
      call device_invcol2(D(2)%x_d, coef%B_d, n)
      call device_invcol2(D(3)%x_d, coef%B_d, n)
    else
      call invcol2(D(1)%x, coef%B, n)
      call invcol2(D(2)%x, coef%B, n)
      call invcol2(D(3)%x, coef%B, n)
    end if
    
    call gs%op(D(1), GS_OP_ADD)
    call gs%op(D(2), GS_OP_ADD)
    call gs%op(D(3), GS_OP_ADD)
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_col2(D(1)%x_d, coef%mult_d, n)
      call device_col2(D(2)%x_d, coef%mult_d, n)
      call device_col2(D(3)%x_d, coef%mult_d, n)
    else
      call col2(D(1)%x, coef%mult, n)
      call col2(D(2)%x, coef%mult, n)
      call col2(D(3)%x, coef%mult, n)
    end if
    call field_copy(residual_viscosity(1)%ptr, D(1))
    call field_copy(residual_viscosity(2)%ptr, D(2))
    call field_copy(residual_viscosity(3)%ptr, D(3))
    call field_absval(residual_viscosity(1)%ptr)
    call field_absval(residual_viscosity(2)%ptr)
    call field_absval(residual_viscosity(3)%ptr)

    call field_cmult(residual_viscosity(1)%ptr, this%c_E)
    call field_cmult(residual_viscosity(2)%ptr, this%c_E)
    call field_cmult(residual_viscosity(3)%ptr, this%c_E)

    u_avg = - glsc2(u%x, coef%B, n) / this%volume_domain
    v_avg = - glsc2(v%x, coef%B, n) / this%volume_domain
    w_avg = - glsc2(w%x, coef%B, n) / this%volume_domain
    call field_cadd2(this%u_var, u, u_avg)
    call field_cadd2(this%v_var, v, v_avg)
    call field_cadd2(this%w_var, w, w_avg)
    call field_absval(this%u_var)
    call field_absval(this%v_var)
    call field_absval(this%w_var)

    scaling_uvw(1) = this%c_E / (dt ** ext_bdf%diffusion_time_order) / &
                     glmax(this%u_var%x, u%dof%size())
    scaling_uvw(2) = this%c_E / (dt ** ext_bdf%diffusion_time_order) / &
                     glmax(this%v_var%x, u%dof%size())
    scaling_uvw(3) = this%c_E / (dt ** ext_bdf%diffusion_time_order) / &
                     glmax(this%w_var%x, u%dof%size())

    end associate

    do i = 1, this%n_scalars

       !! estimate by the difference of the extrapolated and the solved advection
       associate(s => this%s(i)%ptr, ext_bdf => this%ext_bdf, &
                 dt => time%dt, coef => this%coef, D => this%D(i+3), &
                 adv => this%adv, gs => this%coef%gs_h, ta => ta(i+3)%ptr, &
                 u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
                 s_avg => s_avg(i), s_var => this%s_var(i), &
                 scaling_s => scaling_s(i), &
                 residual_viscosity => this%residual_viscosity(i+3)%ptr)

       n = s%dof%size()
       call field_rzero(ta)
       call adv%compute_scalar(u, v, w, s, ta, Xh, coef, n)
       call field_sub2(D, ta, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_invcol2(D%x_d, coef%B_d, n)
       else
          call invcol2(D%x, coef%B, n)
       end if
       
       call gs%op(D, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(D%x_d, coef%mult_d, n)
       else
          call col2(D%x, coef%mult, n)
       end if
       call field_absval(D)

       
       s_avg = - glsc2(s%x, coef%B, n) / this%volume_domain
       call field_cadd2(s_var, s, s_avg)
       call field_absval(s_var)

       scaling_s = this%c_E / (dt ** ext_bdf%diffusion_time_order) / &
                     glmax(s_var%x, u%dof%size())

       ! From now on, we use the D field to store the residual viscosity
       call field_cmult(D, scaling_s)
       call field_col2(D, this%h_np2)
       do j = 1, ext_bdf%diffusion_time_order
          call field_invcol2(D, this%vel_mag)
       end do
       end associate
       !! estimate by ds/dt + ui ds/dxi
      !  associate(s => this%s(i)%ptr, ext_bdf => this%ext_bdf, &
      !            dt => time%dt, coef => this%coef, wa => this%wa(i), &
      !            D => this%D(i), gs => this%coef%gs_h, &
      !            adv => this%adv, &
      !            u => this%u, v => this%v, w => this%w, Xh => this%coef%Xh, &
      !            residual_viscosity => this%residual_viscosity(i)%ptr)

      !  n = s%dof%size()

      !  call field_copy(ta, s)
      !  call field_cmult(ta, ext_bdf%diffusion_coeffs(1)/dt)
      !  if (NEKO_BCKND_DEVICE .eq. 1) then
      !     call device_invcol2(wa%x_d, coef%B_d, n)
      !  else
      !     call invcol2(wa%x, coef%B, n)
      !  end if
      !  call field_sub2(ta, wa)
      !  call field_copy(D, ta)

      !  ! advection part
      !  call field_rzero(ta, n)
      !  call adv%compute_scalar(u, v, w, s, ta, &
      !         Xh, coef, n)
      !  if (NEKO_BCKND_DEVICE .eq. 1) then
      !     call device_invcol2(ta%x_d, coef%B_d, n)
      !  else
      !     call invcol2(ta%x, coef%B, n)
      !  end if
      !  call field_sub2(D, ta, n)
      !  call field_copy(residual_viscosity, D)
      !  call field_absval(residual_viscosity)

      !  ! it should be scaled by f(ext_bdf%diffusion_time_order)
      !  ! preliminary, f could be 0.01937*exp(-5.7363*ext_bdf%diffusion_time_order)
      !  ! Could be determined afterwards
      !  call field_cmult(residual_viscosity, &
      !      this%c_E)

      !  end associate
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine residual_viscosity_compute

end module residual_viscosity
