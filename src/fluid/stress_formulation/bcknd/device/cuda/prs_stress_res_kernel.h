/*
 Copyright (c) 2024, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __FLUID_PRS_STRESS_RES_KERNEL__
#define __FLUID_PRS_STRESS_RES_KERNEL__

template< typename T >
__global__ void prs_stress_res_svv_part1_1_kernel(T * __restrict__ ta1,
                                                  T * __restrict__ ta2,
                                                  T * __restrict__ ta3,
                                                  T * __restrict__ wa1,
                                                  T * __restrict__ wa2,
                                                  T * __restrict__ wa3,
                                                  const T * __restrict__ s11,
                                                  const T * __restrict__ s22,
                                                  const T * __restrict__ s33,
                                                  const T * __restrict__ s12,
                                                  const T * __restrict__ s13,
                                                  const T * __restrict__ s23,
                                                  const T * __restrict__ f_u,
                                                  const T * __restrict__ f_v,
                                                  const T * __restrict__ f_w,
                                                  const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    wa1[i] -= 2.0 * (ta1[i] * s11[i] 
                   + ta2[i] * s12[i]
                   + ta3[i] * s13[i]);
    wa2[i] -= 2.0 * (ta1[i] * s12[i]
                   + ta2[i] * s22[i]
                   + ta3[i] * s23[i]);
    wa3[i] -= 2.0 * (ta1[i] * s13[i] 
                   + ta2[i] * s23[i] 
                   + ta3[i] * s33[i]);
  }

}

template< typename T >
__global__ void prs_stress_res_svv_part1_2_kernel(T * __restrict__ work1,
                                                  T * __restrict__ work2,
                                                  T * __restrict__ work3,
                                                  T * __restrict__ ta1,
                                                  T * __restrict__ ta2,
                                                  T * __restrict__ ta3,
                                                  const T * __restrict__ svv_h1,
                                                  const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    work1[i] *= svv_h1[i];
    work2[i] *= svv_h1[i];
    work3[i] *= svv_h1[i];
    ta1[i] *= svv_h1[i];
    ta2[i] *= svv_h1[i];
    ta3[i] *= svv_h1[i];
  }

}

template< typename T >
__global__ void prs_stress_res_svv_part1_3_kernel(T * __restrict__ wa1,
                                                  T * __restrict__ wa2,
                                                  T * __restrict__ wa3,
                                                  const T * __restrict__ a11,
                                                  const T * __restrict__ a12,
                                                  const T * __restrict__ a13,
                                                  const T * __restrict__ a21,
                                                  const T * __restrict__ a22,
                                                  const T * __restrict__ a23,
                                                  const T * __restrict__ a31,
                                                  const T * __restrict__ a32,
                                                  const T * __restrict__ a33,
                                                  const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    wa1[i] -= 2.0 * (a11[i] + a12[i] + a13[i]);
    wa1[i] -= 2.0 * (a21[i] + a22[i] + a23[i]);
    wa1[i] -= 2.0 * (a31[i] + a32[i] + a33[i]);
  }

}

template< typename T >
__global__ void prs_stress_res_svv_part1_4_kernel(T * __restrict__ ta1,
                                              T * __restrict__ ta2,
                                              T * __restrict__ ta3,
                                              T * __restrict__ wa1,
                                              T * __restrict__ wa2,
                                              T * __restrict__ wa3,
                                              const T * __restrict__ f_u,
                                              const T * __restrict__ f_v,
                                              const T * __restrict__ f_w,
                                              const T * __restrict__ B,
                                              const T * __restrict__ rho,
                                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    ta1[i] = (f_u[i] / rho[i]) - ((wa1[i] / rho[i]) * B[i]);
    ta2[i] = (f_v[i] / rho[i]) - ((wa2[i] / rho[i]) * B[i]);
    ta3[i] = (f_w[i] / rho[i]) - ((wa3[i] / rho[i]) * B[i]);
  }

}


template< typename T >
__global__ void prs_stress_res_part1_kernel(T * __restrict__ ta1,
                                            T * __restrict__ ta2,
                                            T * __restrict__ ta3,
                                            T * __restrict__ wa1,
                                            T * __restrict__ wa2,
                                            T * __restrict__ wa3,
                                            const T * __restrict__ s11,
                                            const T * __restrict__ s22,
                                            const T * __restrict__ s33,
                                            const T * __restrict__ s12,
                                            const T * __restrict__ s13,
                                            const T * __restrict__ s23,
                                            const T * __restrict__ f_u,
                                            const T * __restrict__ f_v,
                                            const T * __restrict__ f_w,
                                            const T * __restrict__ B,
                                            const T * __restrict__ rho,
                                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    wa1[i] -= 2.0 * (ta1[i] * s11[i] 
                   + ta2[i] * s12[i]
                   + ta3[i] * s13[i]);
    wa2[i] -= 2.0 * (ta1[i] * s12[i]
                   + ta2[i] * s22[i]
                   + ta3[i] * s23[i]);
    wa3[i] -= 2.0 * (ta1[i] * s13[i] 
                   + ta2[i] * s23[i] 
                   + ta3[i] * s33[i]);

    ta1[i] = (f_u[i] / rho[i]) - ((wa1[i] / rho[i]) * B[i]);
    ta2[i] = (f_v[i] / rho[i]) - ((wa2[i] / rho[i]) * B[i]);
    ta3[i] = (f_w[i] / rho[i]) - ((wa3[i] / rho[i]) * B[i]);
  }

}

template< typename T >
__global__ void prs_stress_res_part3_kernel(T * __restrict__ p_res,
                                            const T * __restrict__ ta1,
                                            const T * __restrict__ ta2,
                                            const T * __restrict__ ta3,
                                            const T * __restrict__ wa1,
                                            const T * __restrict__ wa2,
                                            const T * __restrict__ wa3,
                                            const T dtbd,
                                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    p_res[i] = p_res[i] - (dtbd * (ta1[i] + ta2[i] + ta3[i]))
      - (wa1[i] + wa2[i] + wa3[i]);
  }
}

#endif // __FLUID_PRS_STRESS_RES_KERNEL__
