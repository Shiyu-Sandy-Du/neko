! Copyright (c) 2024, The Neko Authors
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
!> Implements the cpu kernel for the `strain_rate_based_stress_t` type.
!! Maintainer: Shiyu Du.

module strain_rate_based_stress_cpu
  use num_types, only : rp
  use field_list, only : field_list_t
  use math, only : vcross
  use field, only : field_t
  use operators, only : dudxyz, strain_rate
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use gs_ops, only : GS_OP_ADD
  use coefs, only : coef_t
  implicit none
  private

  public :: strain_rate_based_stress_compute_cpu

contains

  !> Computes the generic Coriolis source term on the cpu.
  !! @param fields The right-hand side, which should be the velocity components.
  !! @param omega The rotation vector.
  !! @param omega The geostrophic wind.
  subroutine strain_rate_based_stress_compute_cpu(fields, nut, coef)
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(in) :: nut
    type(coef_t), intent(in) :: coef
    integer :: i, n
    type(field_t), pointer :: f_x, f_y, f_z, u, v, w
    type(field_t), pointer :: s11, s22, s33, s12, s13, s23
    type(field_t), pointer :: tauijd1, tauijd2, tauijd3
    integer :: temp_indices(9)

    n = fields%item_size(1)

    f_x => fields%get_by_index(1)
    f_y => fields%get_by_index(2)
    f_z => fields%get_by_index(3)

    call neko_scratch_registry%request_field(s11, temp_indices(1))
    call neko_scratch_registry%request_field(s22, temp_indices(2))
    call neko_scratch_registry%request_field(s33, temp_indices(3))
    call neko_scratch_registry%request_field(s12, temp_indices(4))
    call neko_scratch_registry%request_field(s13, temp_indices(5))
    call neko_scratch_registry%request_field(s23, temp_indices(6))
    call neko_scratch_registry%request_field(tauijd1, temp_indices(7))
    call neko_scratch_registry%request_field(tauijd2, temp_indices(8))
    call neko_scratch_registry%request_field(tauijd3, temp_indices(9))

    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("w")

    call strain_rate(s11%x, s22%x, s33%x, s12%x, s13%x, s23%x, u, v, w, coef)

    call coef%gs_h%op(s11, GS_OP_ADD)
    call coef%gs_h%op(s22, GS_OP_ADD)
    call coef%gs_h%op(s33, GS_OP_ADD)
    call coef%gs_h%op(s12, GS_OP_ADD)
    call coef%gs_h%op(s13, GS_OP_ADD)
    call coef%gs_h%op(s23, GS_OP_ADD)

    ! make stress tensor tau_ij from s_ij
    do concurrent (i = 1:n)
       s11%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s11%x(i,1,1,1) * coef%mult(i,1,1,1)
       s22%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s22%x(i,1,1,1) * coef%mult(i,1,1,1)
       s33%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s33%x(i,1,1,1) * coef%mult(i,1,1,1)
       s12%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s12%x(i,1,1,1) * coef%mult(i,1,1,1)
       s13%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s13%x(i,1,1,1) * coef%mult(i,1,1,1)
       s23%x(i,1,1,1) = 2.0_rp * nut%x(i,1,1,1) * s23%x(i,1,1,1) * coef%mult(i,1,1,1)
    end do

    ! d tau_ij / d x_j
    call dudxyz(tauijd1%x, s11%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(tauijd2%x, s12%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(tauijd3%x, s13%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    do i = 1,n
       f_x%x(i,1,1,1) = tauijd1%x(i,1,1,1) + &
                        tauijd2%x(i,1,1,1) + tauijd3%x(i,1,1,1)
    end do

    call dudxyz(tauijd1%x, s12%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(tauijd2%x, s22%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(tauijd3%x, s23%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    do i = 1,n
       f_y%x(i,1,1,1) = tauijd1%x(i,1,1,1) + &
                        tauijd2%x(i,1,1,1) + tauijd3%x(i,1,1,1)
    end do

    call dudxyz(tauijd1%x, s13%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call dudxyz(tauijd2%x, s23%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call dudxyz(tauijd3%x, s33%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    do i = 1,n
       f_z%x(i,1,1,1) = tauijd1%x(i,1,1,1) + &
                        tauijd2%x(i,1,1,1) + tauijd3%x(i,1,1,1)
    end do

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine strain_rate_based_stress_compute_cpu

end module strain_rate_based_stress_cpu
