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
module ax_helm_svv_diffcomp_cpu
  use ax_helm_svv_diffcomp, only : ax_helm_svv_diffcomp_t
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
  type, public, extends(ax_helm_svv_diffcomp_t) :: ax_helm_svv_diffcomp_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute_vector => &
                                 ax_helm_svv_diffcomp_compute_vector
  end type ax_helm_svv_diffcomp_cpu_t

contains

  !> Compute \f$ Ax \f$ inside a Krylov method, taking 3
  !! components of a vector field in a coupled manner.
  !! @param au Result for the first component of the vector.
  !! @param av Result for the first component of the vector.
  !! @param aw Result for the first component of the vector.
  !! @param u The first component of the vector.
  !! @param v The second component of the vector.
  !! @param w The third component of the vector.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  subroutine ax_helm_svv_diffcomp_compute_vector(this, au, av, aw, u, v, w, &
             coef, msh, Xh)
    class(ax_helm_svv_diffcomp_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)


    call ax_helm_svv_diffcomp_lx(au, av, aw, u, v, w, &
            Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, &
            coef%dsdy, coef%dsdz, coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, this%svv%h1_1, this%svv%h1_2, this%svv%h1_3, &
            this%svv%filter%fh, this%svv%filter%fht, msh%nelv, Xh%lx)

    if (coef%ifh2) then
       call addcol4 (au, coef%h2, coef%B, u, coef%dof%size())
       call addcol4 (av, coef%h2, coef%B, v, coef%dof%size())
       call addcol4 (aw, coef%h2, coef%B, w, coef%dof%size())
    end if
  end subroutine ax_helm_svv_diffcomp_compute_vector



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
  subroutine ax_helm_svv_diffcomp_lx(au, av, aw, u, v, w, &
       Dx, Dy, Dz, Dxt, Dyt, Dzt, h1, &
       drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, svv_h1_1, svv_h1_2, svv_h1_3, svv_Q, svv_Qt, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
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
    real(kind=rp), intent(in) :: svv_h1_1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_h1_2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_h1_3(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_Q(lx, lx), svv_Qt(lx, lx)
    real(kind=rp) :: ur_h(lx, lx, lx)
    real(kind=rp) :: us_h(lx, lx, lx)
    real(kind=rp) :: ut_h(lx, lx, lx)
    real(kind=rp) :: vr_h(lx, lx, lx)
    real(kind=rp) :: vs_h(lx, lx, lx)
    real(kind=rp) :: vt_h(lx, lx, lx)
    real(kind=rp) :: wr_h(lx, lx, lx)
    real(kind=rp) :: ws_h(lx, lx, lx)
    real(kind=rp) :: wt_h(lx, lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
   !  real(kind=rp) :: u_svv(lx, lx, lx, n)
    real(kind=rp) :: ur_svv(lx, lx, lx)
    real(kind=rp) :: us_svv(lx, lx, lx)
    real(kind=rp) :: ut_svv(lx, lx, lx)
    real(kind=rp) :: vr_svv(lx, lx, lx)
    real(kind=rp) :: vs_svv(lx, lx, lx)
    real(kind=rp) :: vt_svv(lx, lx, lx)
    real(kind=rp) :: wr_svv(lx, lx, lx)
    real(kind=rp) :: ws_svv(lx, lx, lx)
    real(kind=rp) :: wt_svv(lx, lx, lx)
    
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + Dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + Dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + Dx(i,k) * w(k,j,1,e)
             end do
             wur(i,j,1) = tmp1
             wvr(i,j,1) = tmp2
             wwr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + Dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + Dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + Dy(j,l) * w(i,l,k,e)
                end do
                wus(i,j,k) = tmp1
                wvs(i,j,k) = tmp2
                wws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + Dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + Dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + Dz(k,l) * w(i,1,l,e)
             end do
             wut(i,1,k) = tmp1
             wvt(i,1,k) = tmp2
             wwt(i,1,k) = tmp3
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
          vr(i,1,1) = (drdx(i,1,1,e) * wvr(i,1,1) &
                     + dsdx(i,1,1,e) * wvs(i,1,1) &
                     + dtdx(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          vs(i,1,1) = (drdy(i,1,1,e) * wvr(i,1,1) &
                     + dsdy(i,1,1,e) * wvs(i,1,1) &
                     + dtdy(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          vt(i,1,1) = (drdz(i,1,1,e) * wvr(i,1,1) &
                     + dsdz(i,1,1,e) * wvs(i,1,1) &
                     + dtdz(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          wr(i,1,1) = (drdx(i,1,1,e) * wwr(i,1,1) &
                     + dsdx(i,1,1,e) * wws(i,1,1) &
                     + dtdx(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
          ws(i,1,1) = (drdy(i,1,1,e) * wwr(i,1,1) &
                     + dsdy(i,1,1,e) * wws(i,1,1) &
                     + dtdy(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
          wt(i,1,1) = (drdz(i,1,1,e) * wwr(i,1,1) &
                     + dsdz(i,1,1,e) * wws(i,1,1) &
                     + dtdz(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
       end do

       ! spatial convolution for spectral vanishing
       call tnsr3d_el(ur_svv, lx, ur, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(us_svv, lx, us, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(ut_svv, lx, ut, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(vr_svv, lx, vr, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(vs_svv, lx, vs, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(vt_svv, lx, vt, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(wr_svv, lx, wr, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(ws_svv, lx, ws, lx, svv_Q, svv_Qt, svv_Qt)
       call tnsr3d_el(wt_svv, lx, wt, lx, svv_Q, svv_Qt, svv_Qt)

       do i = 1, lx*lx*lx
          ! multiply the viscosity
          ur_h(i,1,1) = (svv_h1_1(i,1,1,e) * ur_svv(i,1,1) + &
                        h1(i,1,1,e) * ur(i,1,1)) * weights3(i,1,1)
          us_h(i,1,1) = (svv_h1_1(i,1,1,e) * us_svv(i,1,1) + &
                        h1(i,1,1,e) * us(i,1,1)) * weights3(i,1,1)
          ut_h(i,1,1) = (svv_h1_1(i,1,1,e) * ut_svv(i,1,1) + &
                        h1(i,1,1,e) * ut(i,1,1)) * weights3(i,1,1)
          vr_h(i,1,1) = (svv_h1_2(i,1,1,e) * vr_svv(i,1,1) + &
                        h1(i,1,1,e) * vr(i,1,1)) * weights3(i,1,1)
          vs_h(i,1,1) = (svv_h1_2(i,1,1,e) * vs_svv(i,1,1) + &
                        h1(i,1,1,e) * vs(i,1,1)) * weights3(i,1,1)
          vt_h(i,1,1) = (svv_h1_2(i,1,1,e) * vt_svv(i,1,1) + &
                        h1(i,1,1,e) * vt(i,1,1)) * weights3(i,1,1)
          wr_h(i,1,1) = (svv_h1_3(i,1,1,e) * wr_svv(i,1,1) + &
                        h1(i,1,1,e) * wr(i,1,1)) * weights3(i,1,1)
          ws_h(i,1,1) = (svv_h1_3(i,1,1,e) * ws_svv(i,1,1) + &
                        h1(i,1,1,e) * ws(i,1,1)) * weights3(i,1,1)
          wt_h(i,1,1) = (svv_h1_3(i,1,1,e) * wt_svv(i,1,1) + &
                        h1(i,1,1,e) * wt(i,1,1)) * weights3(i,1,1)
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
          wvr(i,1,1) = drdx(i,1,1,e) * vr_h(i,1,1) &
                     + drdy(i,1,1,e) * vs_h(i,1,1) &
                     + drdz(i,1,1,e) * vt_h(i,1,1)
          wvs(i,1,1) = dsdx(i,1,1,e) * vr_h(i,1,1) &
                     + dsdy(i,1,1,e) * vs_h(i,1,1) &
                     + dsdz(i,1,1,e) * vt_h(i,1,1)
          wvt(i,1,1) = dtdx(i,1,1,e) * vr_h(i,1,1) &
                     + dtdy(i,1,1,e) * vs_h(i,1,1) &
                     + dtdz(i,1,1,e) * vt_h(i,1,1)
          wwr(i,1,1) = drdx(i,1,1,e) * wr_h(i,1,1) &
                     + drdy(i,1,1,e) * ws_h(i,1,1) &
                     + drdz(i,1,1,e) * wt_h(i,1,1)
          wws(i,1,1) = dsdx(i,1,1,e) * wr_h(i,1,1) &
                     + dsdy(i,1,1,e) * ws_h(i,1,1) &
                     + dsdz(i,1,1,e) * wt_h(i,1,1)
          wwt(i,1,1) = dtdx(i,1,1,e) * wr_h(i,1,1) &
                     + dtdy(i,1,1,e) * ws_h(i,1,1) &
                     + dtdz(i,1,1,e) * wt_h(i,1,1)
       end do

       do j = 1, lx*lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + Dxt(i,k) * wur(k,j,1)
                tmp2 = tmp2 + Dxt(i,k) * wvr(k,j,1)
                tmp3 = tmp3 + Dxt(i,k) * wwr(k,j,1)
             end do
             au(i,j,1,e) = tmp1
             av(i,j,1,e) = tmp2
             aw(i,j,1,e) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + Dyt(j,l) * wus(i,l,k)
                   tmp2 = tmp2 + Dyt(j,l) * wvs(i,l,k)
                   tmp3 = tmp3 + Dyt(j,l) * wws(i,l,k)
                end do
                au(i,j,k,e) = au(i,j,k,e) + tmp1
                av(i,j,k,e) = av(i,j,k,e) + tmp2
                aw(i,j,k,e) = aw(i,j,k,e) + tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + Dzt(k,l) * wut(i,1,l)
                tmp2 = tmp2 + Dzt(k,l) * wvt(i,1,l)
                tmp3 = tmp3 + Dzt(k,l) * wwt(i,1,l)
             end do
             au(i,1,k,e) = au(i,1,k,e) + tmp1
             av(i,1,k,e) = av(i,1,k,e) + tmp2
             aw(i,1,k,e) = aw(i,1,k,e) + tmp3
          end do
       end do

    end do
  end subroutine ax_helm_svv_diffcomp_lx

end module ax_helm_svv_diffcomp_cpu
