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
module ax_helm_svv_cpu
  use ax_helm_svv, only : ax_helm_svv_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  use tensor, only : tnsr3d_el, tnsr3d
  use spectral_vanishing_viscosity, only : svv_t
  implicit none
  private

  !> CPU matrix-vector product for a Helmholtz problem.
  type, public, extends(ax_helm_svv_t) :: ax_helm_svv_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute => ax_helm_svv_compute
  end type ax_helm_svv_cpu_t

contains

  !> Compute the product.
  !! @param w Vector of size @a (lx,ly,lz,nelv).
  !! @param u Vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  !! @note Since this is a performance-crtical routine, it is implemented in
  !! several kernels corresponding to different polynmial orders.
  subroutine ax_helm_svv_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_svv_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    integer :: i

    call ax_helm_svv_lx(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, &
            coef%dsdz, coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, this%svv%h1, this%svv%filter%fh, &
            this%svv%filter%fht, msh%nelv, Xh%lx)

    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%size())


  end subroutine ax_helm_svv_compute

  !> Generic CPU kernel for the Helmholz matrix-vector product.
  !! @param w Result.
  !! @param u Velocity field.
  !! @param Dx Derivative operator in first dimension.
  !! @param Dy Derivative operator in second dimension.
  !! @param Dz Derivative operator in third dimension.
  !! @param Dxt Derivative operator transpose in first dimension.
  !! @param Dyt Derivative operator transpose in second dimension.
  !! @param Dzt Derivative operator transpose in third dimension.
  !! @param G11 Geometric factor.
  !! @param n Number of elements.
  !! @param lx Polynomial order.
  subroutine ax_helm_svv_lx(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, svv_h1, svv_Q, svv_Qt, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp), intent(in) :: svv_h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_Q(lx, lx), svv_Qt(lx, lx)
    real(kind=rp) :: ur_h(lx, lx, lx)
    real(kind=rp) :: us_h(lx, lx, lx)
    real(kind=rp) :: ut_h(lx, lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
   !  real(kind=rp) :: u_svv(lx, lx, lx, n)
    real(kind=rp) :: ur_svv(lx, lx, lx)
    real(kind=rp) :: us_svv(lx, lx, lx)
    real(kind=rp) :: ut_svv(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: ident(lx, lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l

    do i = 1, lx
       do j = 1, lx
          if (i .eq. j) then
             ident(i,j) = 1.0_rp
          else
             ident(i,j) = 0.0_rp
          end if
       end do
    end do

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dx(i,k) * u(k,j,1,e)
             end do
             wur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dy(j,l) * u(i,l,k,e)
                end do
                wus(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dz(k,l) * u(i,1,l,e)
             end do
             wut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx*lx*lx
          ur(i,1,1) = (drdx(i,1,1,e) * wur(i,1,1) &
                     + dsdx(i,1,1,e) * wus(i,1,1) &
                     + dtdx(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          us(i,1,1) = (drdy(i,1,1,e) * wur(i,1,1) &
                     + dsdy(i,1,1,e) * wus(i,1,1) &
                     + dtdy(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          ut(i,1,1) = (drdz(i,1,1,e) * wur(i,1,1) &
                     + dsdz(i,1,1,e) * wus(i,1,1) &
                     + dtdz(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
       end do

      !  ! spatial convolution for spectral vanishing
       call tnsr3d_el(ur_svv, lx, ur, lx, svv_Q, svv_Qt, ident)
       call tnsr3d_el(us_svv, lx, us, lx, svv_Q, svv_Qt, ident)
       call tnsr3d_el(ut_svv, lx, ut, lx, svv_Q, svv_Qt, ident)

       do i = 1, lx*lx*lx
          ! multiply the viscosity
          ur_h(i,1,1) = (svv_h1(i,1,1,e) * ur_svv(i,1,1) + &
                        h1(i,1,1,e) * ur(i,1,1)) * weights3(i,1,1)
          us_h(i,1,1) = (svv_h1(i,1,1,e) * us_svv(i,1,1) + &
                        h1(i,1,1,e) * us(i,1,1)) * weights3(i,1,1)
          ut_h(i,1,1) = (svv_h1(i,1,1,e) * ut_svv(i,1,1) + &
                        h1(i,1,1,e) * ut(i,1,1)) * weights3(i,1,1)
          ! utilize wur, wus, wut as work arrays again
          wur(i,1,1) = drdx(i,1,1,e) * ur_h(i,1,1) &
                     + drdy(i,1,1,e) * us_h(i,1,1) &
                     + drdz(i,1,1,e) * ut_h(i,1,1)
          wus(i,1,1) = dsdx(i,1,1,e) * ur_h(i,1,1) &
                     + dsdy(i,1,1,e) * us_h(i,1,1) &
                     + dsdz(i,1,1,e) * ut_h(i,1,1)
          wut(i,1,1) = dtdx(i,1,1,e) * ur_h(i,1,1) &
                     + dtdy(i,1,1,e) * us_h(i,1,1) &
                     + dtdz(i,1,1,e) * ut_h(i,1,1)
       end do

       do j = 1, lx*lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dxt(i,k) * wur(k,j,1)
             end do
             w(i,j,1,e) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dyt(j,l) * wus(i,l,k)
                end do
                w(i,j,k,e) = w(i,j,k,e) + tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dzt(k,l) * wut(i,1,l)
             end do
             w(i,1,k,e) = w(i,1,k,e) + tmp
          end do
       end do

    end do
  end subroutine ax_helm_svv_lx

end module ax_helm_svv_cpu
