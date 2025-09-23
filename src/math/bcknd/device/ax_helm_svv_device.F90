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
!   * Neither the name of the wthors nor the names of its
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
module ax_helm_svv_device
  use ax_helm_svv, only : ax_helm_svv_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use device_math, only : device_addcol4
  use device, only : device_get_ptr, device_alloc
  use tensor_device, only : tnsr3d_device
  use num_types, only : rp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_char, &
                                          c_size_t, C_NULL_PTR, &
                                          c_associated
  implicit none
  private

  type, public, extends(ax_helm_svv_t) :: ax_helm_svv_device_t
     type(c_ptr) :: ur_d = C_NULL_PTR
     type(c_ptr) :: us_d = C_NULL_PTR
     type(c_ptr) :: ut_d = C_NULL_PTR
     type(c_ptr) :: ur_svv_d = C_NULL_PTR
     type(c_ptr) :: us_svv_d = C_NULL_PTR
     type(c_ptr) :: ut_svv_d = C_NULL_PTR
   contains
     procedure, pass(this) :: compute => ax_helm_svv_device_compute
  end type ax_helm_svv_device_t

#ifdef HAVE_CUDA
  interface
     subroutine cuda_ax_helm_svv_part1(ur_d, us_d, ut_d,&
          u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, drdy_d, drdz_d, &
          dsdx_d, dsdy_d, dsdz_d, &
          dtdx_d, dtdy_d, dtdz_d, &
          jacinv_d, nelv, lx) &
          bind(c, name='cuda_ax_helm_svv_part1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ur_d, us_d, ut_d
       type(c_ptr), value :: u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, drdy_d, drdz_d
       type(c_ptr), value :: dsdx_d, dsdy_d, dsdz_d
       type(c_ptr), value :: dtdx_d, dtdy_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nelv, lx
     end subroutine cuda_ax_helm_svv_part1
  end interface
  interface
     subroutine cuda_ax_helm_svv_part2(w_d, &
          ur_d, us_d, ut_d,&
          ur_svv_d, us_svv_d, ut_svv_d, &
          dx_d, dy_d, dz_d, &
          h1_d, drdx_d, drdy_d, drdz_d, &
          dsdx_d, dsdy_d, dsdz_d, &
          dtdx_d, dtdy_d, dtdz_d, &
          w3_d, svv_h1_d, nelv, lx) &
          bind(c, name='cuda_ax_helm_svv_part2')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: w_d
       type(c_ptr), value :: ur_d, us_d, ut_d
       type(c_ptr), value :: ur_svv_d, us_svv_d, ut_svv_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: h1_d
       type(c_ptr), value :: drdx_d, drdy_d, drdz_d
       type(c_ptr), value :: dsdx_d, dsdy_d, dsdz_d
       type(c_ptr), value :: dtdx_d, dtdy_d, dtdz_d
       type(c_ptr), value :: w3_d
       type(c_ptr), value :: svv_h1_d
       integer(c_int) :: nelv, lx
     end subroutine cuda_ax_helm_svv_part2
  end interface
#endif

contains

  subroutine ax_helm_svv_device_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_svv_device_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(c_ptr) :: u_d, w_d
    integer :: n

    associate(lx => Xh%lx, ly => Xh%ly, lz => Xh%lz, &
          nelv => msh%nelv, &
          ur_d => this%ur_d, us_d => this%us_d, ut_d => this%ut_d, &
          ur_svv_d => this%ur_svv_d, us_svv_d => this%us_svv_d, &
          ut_svv_d => this%ut_svv_d, &
          svv_Q => this%svv%filter%fh_d, svv_Qt => this%svv%filter%fht_d, &
          svv_direction => this%svv%direction, &
          ident => this%svv%filter%ident_d)

    u_d = device_get_ptr(u)
    w_d = device_get_ptr(w)
    
#ifdef HAVE_HIP
    call neko_error('HIP is not implemented for SVV')
#elif HAVE_CUDA
    call cuda_ax_helm_svv_part1(ur_d, us_d, ut_d, &
          u_d, &
          Xh%dx_d, Xh%dy_d, Xh%dz_d, &
          coef%drdx_d, coef%drdy_d, coef%drdz_d, &
          coef%dsdx_d, coef%dsdy_d, coef%dsdz_d, &
          coef%dtdx_d, coef%dtdy_d, coef%dtdz_d, &
          coef%jacinv_d, &
          nelv, lx)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for SVV')
#endif

    ! Filtering operation
    if (svv_direction .eq. "rst") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, svv_Q, svv_Qt, svv_Qt, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, svv_Q, svv_Qt, svv_Qt, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, svv_Q, svv_Qt, svv_Qt, nelv)
    else if (svv_direction .eq. "rs") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, svv_Q, svv_Qt, ident, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, svv_Q, svv_Qt, ident, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, svv_Q, svv_Qt, ident, nelv)
    else if (svv_direction .eq. "rt") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, svv_Q, ident, svv_Qt, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, svv_Q, ident, svv_Qt, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, svv_Q, ident, svv_Qt, nelv)
    else if (svv_direction .eq. "st") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, ident, svv_Qt, svv_Qt, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, ident, svv_Qt, svv_Qt, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, ident, svv_Qt, svv_Qt, nelv)
    else if (svv_direction .eq. "r") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, svv_Q, ident, ident, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, svv_Q, ident, ident, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, svv_Q, ident, ident, nelv)
    else if (svv_direction .eq. "s") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, ident, svv_Qt, ident, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, ident, svv_Qt, ident, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, ident, svv_Qt, ident, nelv)
    else if (svv_direction .eq. "t") then
       call tnsr3d_device(ur_svv_d, lx, ur_d, lx, ident, ident, svv_Qt, nelv)
       call tnsr3d_device(us_svv_d, lx, us_d, lx, ident, ident, svv_Qt, nelv)
       call tnsr3d_device(ut_svv_d, lx, ut_d, lx, ident, ident, svv_Qt, nelv)
    end if

#ifdef HAVE_HIP
    call neko_error('HIP is not implemented for SVV')
#elif HAVE_CUDA
    call cuda_ax_helm_svv_part2(w_d, &
          ur_d, us_d, ut_d, &
          ur_svv_d, us_svv_d, ut_svv_d, &
          Xh%dx_d, Xh%dy_d, Xh%dz_d, &
          coef%h1_d, coef%drdx_d, coef%drdy_d, coef%drdz_d, &
          coef%dsdx_d, coef%dsdy_d, coef%dsdz_d, &
          coef%dtdx_d, coef%dtdy_d, coef%dtdz_d, &
          Xh%w3_d, this%svv%h1_d, msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call neko_error('OPENCL is not implemented for SVV')
#endif

    if (coef%ifh2) then
       call device_addcol4(w_d ,coef%h2_d, coef%B_d, u_d, coef%dof%size())
    end if

    end associate

  end subroutine ax_helm_svv_device_compute
end module ax_helm_svv_device
