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
!> Subroutines to add advection terms to the RHS of a transport equation.
module adv_full_dealias
  use advection, only: advection_t
  use num_types, only: rp
  use math, only: vdot3, sub2
  use space, only: space_t, GL
  use field, only: field_t
  use gather_scatter, only: gs_t
  use gs_ops, only : GS_OP_ADD
  use dofmap , only : dofmap_t
  use coefs, only: coef_t
  use device_math, only: device_vdot3, device_sub2
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM, &
    NEKO_BCKND_OPENCL, NEKO_BCKND_CUDA, NEKO_BCKND_HIP
  use operators, only: opgrad, dudxyz
  use interpolation, only: interpolator_t
  use device, only: device_map, device_get_ptr
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  use utils, only: neko_error
  use filter, only: filter_t
  use field_math, only : field_copy
  implicit none
  private

  !> Type encapsulating advection routines with dealiasing
  type, public, extends(advection_t) :: adv_full_dealias_t
     !> Coeffs of the higher-order space
     type(coef_t) :: coef_GL
     !> gs of the higher-order space
     type(gs_t) :: gs_GL
     !> Dofmap of the higher-order space
     type(dofmap_t) :: dm_GL
     !> Coeffs of the original space in the simulation
     type(coef_t), pointer :: coef_GLL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t) :: GLL_to_GL
     !> The additional higher-order space used in dealiasing
     type(space_t) :: Xh_GL
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     real(kind=rp), allocatable :: temp(:), tbf(:)
     !> Temporary arrays
     real(kind=rp), allocatable :: vr(:), vs(:), vt(:)
     !> Non-linear terms
     ! the following array should be global
     ! and stored as field_t for filters
     type(field_t) :: tx, ty, tz, wa
     type(field_t) :: t11, t22, t33, t12, t13, t23
     type(field_t) :: dump_field, dump_field2
     type(field_t) :: tx11, ty22, tz33, &
                      tx12, ty12, tx13, &
                      tz13, ty23, tz23
     !> Device pointer for `temp`
     type(c_ptr) :: temp_d = C_NULL_PTR
     !> Device pointer for `tbf`
     type(c_ptr) :: tbf_d = C_NULL_PTR
     !> Device pointer for `vr`
     type(c_ptr) :: vr_d = C_NULL_PTR
     !> Device pointer for `vs`
     type(c_ptr) :: vs_d = C_NULL_PTR
     !> Device pointer for `vt`
     type(c_ptr) :: vt_d = C_NULL_PTR

   contains
     !> Add the advection term for the fluid, i.e. \f$u \cdot \nabla u \f$, to
     !! the RHS.
     procedure, pass(this) :: compute => compute_advection_full_dealias
     !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to
     !! the RHS.
     procedure, pass(this) :: compute_scalar => compute_scalar_advection_full_dealias
     !> Constructor
     procedure, pass(this) :: init => init_full_dealias
     !> Destructor
     procedure, pass(this) :: free => free_full_dealias
  end type adv_full_dealias_t

contains

  !> Constructor
  !! @param lxd The polynomial order of the space used in the dealiasing.
  !! @param coef The coefficients of the (space, mesh) pair.
  subroutine init_full_dealias(this, lxd, coef)
    class(adv_full_dealias_t), target, intent(inout) :: this
    integer, intent(in) :: lxd
    type(coef_t), intent(inout), target :: coef
    integer :: nel, n_GL, n

    call this%Xh_GL%init(GL, lxd, lxd, lxd)
    this%Xh_GLL => coef%Xh
    this%coef_GLL => coef
    call this%GLL_to_GL%init(this%Xh_GL, this%Xh_GLL)
    
    call this%dm_GL%init(coef%msh, this%Xh_GL)
    call this%gs_GL%init(this%dm_GL)
    call this%coef_GL%init(this%gs_GL)

    nel = coef%msh%nelv
    n_GL = nel*this%Xh_GL%lxyz
    n = nel*coef%Xh%lxyz
    call this%GLL_to_GL%map(this%coef_GL%drdx, coef%drdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdx, coef%dsdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdx, coef%dtdx, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdy, coef%drdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdy, coef%dsdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdy, coef%dtdy, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%drdz, coef%drdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dsdz, coef%dsdz, nel, this%Xh_GL)
    call this%GLL_to_GL%map(this%coef_GL%dtdz, coef%dtdz, nel, this%Xh_GL)
    if ((NEKO_BCKND_HIP .eq. 1) .or. (NEKO_BCKND_CUDA .eq. 1) .or. &
       (NEKO_BCKND_OPENCL .eq. 1) .or. (NEKO_BCKND_SX .eq. 1) .or. &
       (NEKO_BCKND_XSMM .eq. 1)) then
       allocate(this%temp(n_GL))
       allocate(this%tbf(n_GL))
       allocate(this%vr(n_GL))
       allocate(this%vs(n_GL))
       allocate(this%vt(n_GL))
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%temp, this%temp_d, n_GL)
       call device_map(this%tbf, this%tbf_d, n_GL)
       call device_map(this%vr, this%vr_d, n_GL)
       call device_map(this%vs, this%vs_d, n_GL)
       call device_map(this%vt, this%vt_d, n_GL)
    end if

    call this%tx%init(this%dm_GL, "tx")
    call this%ty%init(this%dm_GL, "ty")
    call this%tz%init(this%dm_GL, "tz")
    call this%wa%init(this%dm_GL, "wa")
    call this%t11%init(this%dm_GL, "t11")
    call this%t22%init(this%dm_GL, "t22")
    call this%t33%init(this%dm_GL, "t33")
    call this%t12%init(this%dm_GL, "t12")
    call this%t13%init(this%dm_GL, "t13")
    call this%t23%init(this%dm_GL, "t23")
    call this%dump_field%init(this%dm_GL, "dump_field")
    call this%dump_field2%init(this%dm_GL, "dump_field2")
    call this%tx11%init(this%dm_GL, "tx11")
    call this%ty22%init(this%dm_GL, "ty22")
    call this%tz33%init(this%dm_GL, "tz33")
    call this%tx12%init(this%dm_GL, "tx12")
    call this%ty12%init(this%dm_GL, "ty12")
    call this%tx13%init(this%dm_GL, "tx13")
    call this%tz13%init(this%dm_GL, "tz13")
    call this%ty23%init(this%dm_GL, "ty23")
    call this%tz23%init(this%dm_GL, "tz23")

  end subroutine init_full_dealias

  !> Destructor
  subroutine free_full_dealias(this)
    class(adv_full_dealias_t), intent(inout) :: this
  end subroutine free_full_dealias


  !> Add the advection term for the fluid, i.e. \f$ \frac{\partial u_i u_j}{x_j} \f$, to
  !! the RHS.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param fx The x component of source term.
  !! @param fy The y component of source term.
  !! @param fz The z component of source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_advection_full_dealias(this, vx, vy, vz, fx, fy, fz, Xh, &
                                       coef, n, dt)
    class(adv_full_dealias_t), intent(inout) :: this
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: fx, fy, fz
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: tfx, tfy, tfz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vr, vs, vt
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: tempx, tempy, tempz
    integer :: e, i, idx, nel, n_GL

    nel = coef%msh%nelv
    n_GL = nel * this%Xh_GL%lxyz

    !This is extremely primitive and unoptimized  on the device //Karp
    associate(c_GL => this%coef_GL, tx => this%tx, ty => this%ty, &
              tz => this%tz, t11 => this%t11, t22 => this%t22, &
              t33 => this%t33, t12 => this%t12, t13 => this%t13, &
              t23 => this%t23, tx11 => this%tx11, ty22 => this%ty22, &
              tz33 => this%tz33, &
              tx12 => this%tx12, ty12 => this%ty12, &
              tx13 => this%tx13, tz13 => this%tz13, &
              ty23 => this%ty23, tz23 => this%tz23, &
              dump_field => this%dump_field, dump_field2 => this%dump_field2)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call neko_error("adv_full_dealiasing not implemented for devices")

      else if ((NEKO_BCKND_SX .eq. 1) .or. (NEKO_BCKND_XSMM .eq. 1)) then

         call neko_error("adv_full_dealiasing not implemented for devices")

      else

         call this%GLL_to_GL%map(tx%x, vx%x, coef%msh%nelv, this%Xh_GL)
         call this%GLL_to_GL%map(ty%x, vy%x, coef%msh%nelv, this%Xh_GL)
         call this%GLL_to_GL%map(tz%x, vz%x, coef%msh%nelv, this%Xh_GL)
         do concurrent (i = 1 : t11%dof%size())
            t11%x(i, 1, 1, 1) = tx%x(i, 1, 1, 1) * tx%x(i, 1, 1, 1)
            t22%x(i, 1, 1, 1) = ty%x(i, 1, 1, 1) * ty%x(i, 1, 1, 1)
            t33%x(i, 1, 1, 1) = tz%x(i, 1, 1, 1) * tz%x(i, 1, 1, 1)
            t12%x(i, 1, 1, 1) = tx%x(i, 1, 1, 1) * ty%x(i, 1, 1, 1)
            t13%x(i, 1, 1, 1) = tx%x(i, 1, 1, 1) * tz%x(i, 1, 1, 1)
            t23%x(i, 1, 1, 1) = ty%x(i, 1, 1, 1) * tz%x(i, 1, 1, 1)
         end do

         call field_copy(this%wa, this%t11)
         call this%explicit_filter%apply(this%t11, this%wa)
         call field_copy(this%wa, this%t22)
         call this%explicit_filter%apply(this%t22, this%wa)
         call field_copy(this%wa, this%t33)
         call this%explicit_filter%apply(this%t33, this%wa)
         call field_copy(this%wa, this%t12)
         call this%explicit_filter%apply(this%t12, this%wa)
         call field_copy(this%wa, this%t13)
         call this%explicit_filter%apply(this%t13, this%wa)
         call field_copy(this%wa, this%t23)
         call this%explicit_filter%apply(this%t23, this%wa)
         
         call opgrad(tx11%x, dump_field%x, dump_field2%x, t11%x, c_GL)
         call opgrad(dump_field%x, ty22%x, dump_field2%x, t22%x, c_GL)
         call opgrad(dump_field%x, dump_field2%x, tz33%x, t33%x, c_GL)
         call opgrad(tx12%x, ty12%x, dump_field%x, t12%x, c_GL)
         call opgrad(tx13%x, dump_field%x, tz13%x, t13%x, c_GL)
         call opgrad(dump_field%x, ty23%x, tz23%x, t23%x, c_GL)


         do e = 1, coef%msh%nelv
            do concurrent (i = 1 : this%Xh_GL%lxyz)
               tfx(i) = tx11%x(i, 1, 1, e) + &
                        ty12%x(i, 1, 1, e) + &
                        tz13%x(i, 1, 1, e)
               tfy(i) = tx12%x(i, 1, 1, e) + &
                        ty22%x(i, 1, 1, e) + &
                        tz23%x(i, 1, 1, e)
               tfz(i) = tx13%x(i, 1, 1, e) + &
                        ty23%x(i, 1, 1, e) + &
                        tz33%x(i, 1, 1, e)
            end do

            call this%GLL_to_GL%map(tempx, tfx, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempy, tfy, 1, this%Xh_GLL)
            call this%GLL_to_GL%map(tempz, tfz, 1, this%Xh_GLL)

            idx = (e-1)*this%Xh_GLL%lxyz+1
            do concurrent (i = 0:this%Xh_GLL%lxyz-1)
               fx%x(i+idx,1,1,1) = fx%x(i+idx,1,1,1) - tempx(i+1)
               fy%x(i+idx,1,1,1) = fy%x(i+idx,1,1,1) - tempy(i+1)
               fz%x(i+idx,1,1,1) = fz%x(i+idx,1,1,1) - tempz(i+1)
            end do
         end do
      end if
    end associate

  end subroutine compute_advection_full_dealias

  !> Add the advection term for a scalar, i.e. \f$u \cdot \nabla s \f$, to the
  !! RHS.
  !! @param this The object.
  !! @param vx The x component of velocity.
  !! @param vy The y component of velocity.
  !! @param vz The z component of velocity.
  !! @param s The scalar.
  !! @param fs The source term.
  !! @param Xh The function space.
  !! @param coef The coefficients of the (Xh, mesh) pair.
  !! @param n Typically the size of the mesh.
  !! @param dt Current time-step, not required for this method.
  subroutine compute_scalar_advection_full_dealias(this, vx, vy, vz, s, fs, Xh, &
                                              coef, n, dt)
    class(adv_full_dealias_t), intent(inout) :: this
    type(field_t), intent(inout) :: vx, vy, vz
    type(field_t), intent(inout) :: s
    type(field_t), intent(inout) :: fs
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n
    real(kind=rp), intent(in), optional :: dt

    real(kind=rp), dimension(this%Xh_GL%lxyz) :: vx_GL, vy_GL, vz_GL, s_GL
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: dsdx, dsdy, dsdz
    real(kind=rp), dimension(this%Xh_GL%lxyz) :: f_GL
    integer :: e, i, idx, nel, n_GL
    real(kind=rp), dimension(this%Xh_GLL%lxyz) :: temp

    call neko_error("compute_scalar_advection_full_dealias not implemented")

  end subroutine compute_scalar_advection_full_dealias

end module adv_full_dealias
