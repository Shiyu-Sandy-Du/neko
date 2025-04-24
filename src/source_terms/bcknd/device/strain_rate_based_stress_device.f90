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
!> Implements the device kernel for the `strain_rate_based_stress_t` type.
!! Maintainer: Shiyu Du.

module strain_rate_based_stress_device
  use num_types, only : rp
  use field_list, only : field_list_t
  use device_math, only : device_invcol2, device_cmult2
  use field, only : field_t
  use operators, only : dudxyz, strain_rate
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use gs_ops, only : GS_OP_ADD
  use coefs, only : coef_t
  use ax_product, only: ax_t
  implicit none
  private

  public :: strain_rate_based_stress_compute_device

contains

  !> Computes the generic Coriolis source term on the device.
  !! @param fields The right-hand side, which should be the velocity components.
  !! @param omega The rotation vector.
  !! @param omega The geostrophic wind.
  subroutine strain_rate_based_stress_compute_device(Ax, fields, nut, coef)
    class(ax_t), intent(in) :: Ax
    type(field_list_t), intent(inout) :: fields
    type(field_t), intent(in) :: nut
    type(coef_t), intent(inout) :: coef
    integer :: n
    type(field_t), pointer :: u, v, w
    
    n = fields%item_size(1)

    u => neko_field_registry%get_field_by_name("u")
    v => neko_field_registry%get_field_by_name("v")
    w => neko_field_registry%get_field_by_name("w")

    ! set up Helmholtz operators for the Laplacian
    call device_cmult2(coef%h1_d, nut%x_d, -1.0_rp, n)
    coef%ifh2 = .false.

    call Ax%compute_vector(fields%x(1), fields%x(2), fields%x(3), &
         u%x, v%x, w%x, coef, coef%msh, coef%Xh)

    call device_invcol2(fields%x_d(1), coef%B_d, n)
    call device_invcol2(fields%x_d(2), coef%B_d, n)
    call device_invcol2(fields%x_d(3), coef%B_d, n)

  end subroutine strain_rate_based_stress_compute_device

end module strain_rate_based_stress_device
