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
module ax_helm_diffcomp_cpu
  use ax_helm_diffcomp, only : ax_helm_diffcomp_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  implicit none
  private

  !> CPU matrix-vector product for a Helmholtz problem.
  type, public, extends(ax_helm_diffcomp_t) :: ax_helm_diffcomp_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute_vector => ax_helm_diffcomp_compute_vector
  end type ax_helm_diffcomp_cpu_t

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
  subroutine ax_helm_diffcomp_compute_vector(this, au, av, aw, u, v, w, coef, msh,&
                                         Xh)
    class(ax_helm_diffcomp_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)


    call ax_helm_diffcomp_lx(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1_1, coef%h1_2, coef%h1_3, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv, Xh%lx)

    if (coef%ifh2) then
       call addcol4 (au, coef%h2, coef%B, u, coef%dof%size())
       call addcol4 (av, coef%h2, coef%B, v, coef%dof%size())
       call addcol4 (aw, coef%h2, coef%B, w, coef%dof%size())
    end if


  end subroutine ax_helm_diffcomp_compute_vector

  subroutine ax_helm_diffcomp_lx(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1_1, h1_2, h1_3, G11, G22, G33, G12, G13, G23, n, lx)

    integer, intent(in) :: n, lx
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1_1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1_2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1_3(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)
    real(kind=rp) :: tmp, t1, t2, t3
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dx(i,k) * u(k,j,1,e)
                t2 = t2 + Dx(i,k) * v(k,j,1,e)
                t3 = t3 + Dx(i,k) * w(k,j,1,e)
             end do
             wur(i,j,1) = t1
             wvr(i,j,1) = t2
             wwr(i,j,1) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_rp
                t2 = 0.0_rp
                t3 = 0.0_rp
                do l = 1, lx
                   t1 = t1 + Dy(j,l) * u(i,l,k,e)
                   t2 = t2 + Dy(j,l) * v(i,l,k,e)
                   t3 = t3 + Dy(j,l) * w(i,l,k,e)
                end do
                wus(i,j,k) = t1
                wvs(i,j,k) = t2
                wws(i,j,k) = t3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do l = 1, lx
                t1 = t1 + Dz(k,l) * u(i,1,l,e)
                t2 = t2 + Dz(k,l) * v(i,1,l,e)
                t3 = t3 + Dz(k,l) * w(i,1,l,e)
             end do
             wut(i,1,k) = t1
             wvt(i,1,k) = t2
             wwt(i,1,k) = t3
          end do
       end do

       do i = 1, lx*lx*lx
          ur(i,1,1) = h1_1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1_1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1_1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )
          vr(i,1,1) = h1_2(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1_2(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1_2(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )
          wr(i,1,1) = h1_3(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1_3(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1_3(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dxt(i,k) * ur(k,j,1)
                t2 = t2 + Dxt(i,k) * vr(k,j,1)
                t3 = t3 + Dxt(i,k) * wr(k,j,1)
             end do
             au(i,j,1,e) = t1
             av(i,j,1,e) = t2
             aw(i,j,1,e) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_rp
                t2 = 0.0_rp
                t3 = 0.0_rp
                do l = 1, lx
                   t1 = t1 + Dyt(j,l) * us(i,l,k)
                   t2 = t2 + Dyt(j,l) * vs(i,l,k)
                   t3 = t3 + Dyt(j,l) * ws(i,l,k)
                end do
                au(i,j,k,e) = au(i,j,k,e) + t1
                av(i,j,k,e) = av(i,j,k,e) + t2
                aw(i,j,k,e) = aw(i,j,k,e) + t3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do l = 1, lx
                t1 = t1 + Dzt(k,l) * ut(i,1,l)
                t2 = t2 + Dzt(k,l) * vt(i,1,l)
                t3 = t3 + Dzt(k,l) * wt(i,1,l)
             end do
             au(i,1,k,e) = au(i,1,k,e) + t1
             av(i,1,k,e) = av(i,1,k,e) + t2
             aw(i,1,k,e) = aw(i,1,k,e) + t3
          end do
       end do

    end do
  end subroutine ax_helm_diffcomp_lx

end module ax_helm_diffcomp_cpu
