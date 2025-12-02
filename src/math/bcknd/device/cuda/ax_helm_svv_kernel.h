#ifndef __MATH_AX_HELM_SVV_KERNEL_H__
#define __MATH_AX_HELM_SVV_KERNEL_H__
/*
 Copyright (c) 2025, The Neko Authors
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

/**
 * Device kernels for Ax helm svv
 */

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_part1_kernel_kstep(T * __restrict__ ur,
                                 T * __restrict__ us,
                                 T * __restrict__ ut,
                                 const T * __restrict__ u,
                                 const T * __restrict__ dx,
                                 const T * __restrict__ dy,
                                 const T * __restrict__ dz,
                                 const T * __restrict__ drdx,
                                 const T * __restrict__ drdy,
                                 const T * __restrict__ drdz,
                                 const T * __restrict__ dsdx,
                                 const T * __restrict__ dsdy,
                                 const T * __restrict__ dsdz,
                                 const T * __restrict__ dtdx,
                                 const T * __restrict__ dtdy,
                                 const T * __restrict__ dtdz,
                                 const T * __restrict__ jacinv) {

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  __shared__ T shu[LX * LX];
  
  T ru[LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ele = e*LX*LX*LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
  }

  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T drdx_local = drdx[ijk+ele];
    const T drdy_local = drdy[ijk+ele];
    const T drdz_local = drdz[ijk+ele];
    const T dsdx_local = dsdx[ijk+ele];
    const T dsdy_local = dsdy[ijk+ele];
    const T dsdz_local = dsdz[ijk+ele];
    const T dtdx_local = dtdx[ijk+ele];
    const T dtdy_local = dtdy[ijk+ele];
    const T dtdz_local = dtdz[ijk+ele];
    const T dj = jacinv[ijk+ele];

    T uttmp = 0.0;
    shu[ij] = ru[k];

    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*LX] * ru[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*LX] * shu[l+j*LX];
      ustmp += shdy[j+l*LX] * shu[i+l*LX];
    }
    __syncthreads();
    ur[ijk + ele] = dj * (urtmp * drdx_local + 
                          ustmp * dsdx_local + 
                          uttmp * dtdx_local);
    us[ijk + ele] = dj * (urtmp * drdy_local + 
                          ustmp * dsdy_local + 
                          uttmp * dtdy_local);
    ut[ijk + ele] = dj * (urtmp * drdz_local + 
                          ustmp * dsdz_local + 
                          uttmp * dtdz_local);    
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_part1_kernel_kstep_padded(T * __restrict__ ur,
                                        T * __restrict__ us,
                                        T * __restrict__ ut,
                                        const T * __restrict__ u,
                                        const T * __restrict__ dx,
                                        const T * __restrict__ dy,
                                        const T * __restrict__ dz,
                                        const T * __restrict__ drdx,
                                        const T * __restrict__ drdy,
                                        const T * __restrict__ drdz,
                                        const T * __restrict__ dsdx,
                                        const T * __restrict__ dsdy,
                                        const T * __restrict__ dsdz,
                                        const T * __restrict__ dtdx,
                                        const T * __restrict__ dtdy,
                                        const T * __restrict__ dtdz,
                                        const T * __restrict__ jacinv) {

  __shared__ T shdx[LX * (LX+1)];
  __shared__ T shdy[LX * (LX+1)];
  __shared__ T shdz[LX * (LX+1)];

  __shared__ T shu[LX * (LX+1)];

  T ru[LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*(LX+1);
  const int ele = e*LX*LX*LX;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T drdx_local = drdx[ijk+ele];
    const T drdy_local = drdy[ijk+ele];
    const T drdz_local = drdz[ijk+ele];
    const T dsdx_local = dsdx[ijk+ele];
    const T dsdy_local = dsdy[ijk+ele];
    const T dsdz_local = dsdz[ijk+ele];
    const T dtdx_local = dtdx[ijk+ele];
    const T dtdy_local = dtdy[ijk+ele];
    const T dtdz_local = dtdz[ijk+ele];
    const T dj  = jacinv[ijk+ele];

    T uttmp = 0.0;
    shu[ij_p] = ru[k];

    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*(LX+1)] * ru[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*(LX+1)] * shu[l+j*(LX+1)];
      ustmp += shdy[j+l*(LX+1)] * shu[i+l*(LX+1)];
    }
    __syncthreads();
    ur[ijk + ele] = dj * (urtmp * drdx_local + 
                          ustmp * dsdx_local + 
                          uttmp * dtdx_local);
    us[ijk + ele] = dj * (urtmp * drdy_local + 
                          ustmp * dsdy_local + 
                          uttmp * dtdy_local);
    ut[ijk + ele] = dj * (urtmp * drdz_local + 
                          ustmp * dsdz_local + 
                          uttmp * dtdz_local);
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_part2_kernel_kstep(T * __restrict__ au,
                                 const T * __restrict__ ur,
                                 const T * __restrict__ us,
                                 const T * __restrict__ ut,
                                 const T * __restrict__ ur_svv,
                                 const T * __restrict__ us_svv,
                                 const T * __restrict__ ut_svv,
                                 const T * __restrict__ dx,
                                 const T * __restrict__ dy,
                                 const T * __restrict__ dz,
                                 const T * __restrict__ h1,
                                 const T * __restrict__ drdx,
                                 const T * __restrict__ drdy,
                                 const T * __restrict__ drdz,
                                 const T * __restrict__ dsdx,
                                 const T * __restrict__ dsdy,
                                 const T * __restrict__ dsdz,
                                 const T * __restrict__ dtdx,
                                 const T * __restrict__ dtdy,
                                 const T * __restrict__ dtdz,
                                 const T * __restrict__ w3,
                                 const T * __restrict__ h1_svv) {

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  __shared__ T shur2[LX * LX];
  __shared__ T shus2[LX * LX];
  T rut2;

  T ruw[LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ele = e*LX*LX*LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ruw[k] = 0.0;
  }

  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T drdx_local = drdx[ijk+ele];
    const T drdy_local = drdy[ijk+ele];
    const T drdz_local = drdz[ijk+ele];
    const T dsdx_local = dsdx[ijk+ele];
    const T dsdy_local = dsdy[ijk+ele];
    const T dsdz_local = dsdz[ijk+ele];
    const T dtdx_local = dtdx[ijk+ele];
    const T dtdy_local = dtdy[ijk+ele];
    const T dtdz_local = dtdz[ijk+ele];
    const T dj = w3[ijk]*h1[ijk+ele];
    const T dj_svv = w3[ijk]*h1_svv[ijk+ele];

    T rur = ur[ijk + ele];
    T rus = us[ijk + ele];
    T rut = ut[ijk + ele];
    T rur_svv = rur - ur_svv[ijk + ele];
    T rus_svv = rus - us_svv[ijk + ele];
    T rut_svv = rut - ut_svv[ijk + ele];

    T ur_h = dj * rur + dj_svv * rur_svv;
    T us_h = dj * rus + dj_svv * rus_svv;
    T ut_h = dj * rut + dj_svv * rut_svv;

    shur2[ij] = drdx_local * ur_h +
                drdy_local * us_h +
                drdz_local * ut_h;
    shus2[ij] = dsdx_local * ur_h +
                dsdy_local * us_h +
                dsdz_local * ut_h;
    rut2 =      dtdx_local * ur_h +
                dtdy_local * us_h +
                dtdz_local * ut_h;

    __syncthreads();

    T uwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur2[l+j*LX] * shdx[l+i*LX];
      ruw[l] += rut2 * shdz[k+l*LX];
      uwijke += shus2[i+l*LX] * shdy[l + j*LX];
    }
    __syncthreads();
    ruw[k] += uwijke; 
  }
#pragma unroll
  for (int k = 0; k < LX; ++k){
    au[ij + k*LX*LX + ele] = ruw[k];
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_part2_kernel_kstep_padded(T * __restrict__ au,
                                        const T * __restrict__ ur,
                                        const T * __restrict__ us,
                                        const T * __restrict__ ut,
                                        const T * __restrict__ ur_svv,
                                        const T * __restrict__ us_svv,
                                        const T * __restrict__ ut_svv,
                                        const T * __restrict__ dx,
                                        const T * __restrict__ dy,
                                        const T * __restrict__ dz,
                                        const T * __restrict__ h1,
                                        const T * __restrict__ drdx,
                                        const T * __restrict__ drdy,
                                        const T * __restrict__ drdz,
                                        const T * __restrict__ dsdx,
                                        const T * __restrict__ dsdy,
                                        const T * __restrict__ dsdz,
                                        const T * __restrict__ dtdx,
                                        const T * __restrict__ dtdy,
                                        const T * __restrict__ dtdz,
                                        const T * __restrict__ w3,
                                        const T * __restrict__ h1_svv) {

  __shared__ T shdx[LX * (LX+1)];
  __shared__ T shdy[LX * (LX+1)];
  __shared__ T shdz[LX * (LX+1)];

  __shared__ T shur2[LX * LX];
  __shared__ T shus2[LX * (LX+1)];
  T rut2;
  
  T rur[LX];
  T rus[LX];
  T rut[LX];
  T rur_svv[LX];
  T rus_svv[LX];
  T rut_svv[LX];

  T ruw[LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*(LX+1);
  const int ele = e*LX*LX*LX;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ruw[k] = 0.0;
  }

  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T drdx_local = drdx[ijk+ele];
    const T drdy_local = drdy[ijk+ele];
    const T drdz_local = drdz[ijk+ele];
    const T dsdx_local = dsdx[ijk+ele];
    const T dsdy_local = dsdy[ijk+ele];
    const T dsdz_local = dsdz[ijk+ele];
    const T dtdx_local = dtdx[ijk+ele];
    const T dtdy_local = dtdy[ijk+ele];
    const T dtdz_local = dtdz[ijk+ele];
    const T dj = w3[ijk]*h1[ijk+ele];
    const T dj_svv = w3[ijk]*h1_svv[ijk+ele];

    T rur = ur[ijk + ele];
    T rus = us[ijk + ele];
    T rut = ut[ijk + ele];
    T rur_svv = rur - ur_svv[ijk + ele];
    T rus_svv = rus - us_svv[ijk + ele];
    T rut_svv = rut - ut_svv[ijk + ele];

    T ur_h = dj * rur + dj_svv * rur_svv;
    T us_h = dj * rus + dj_svv * rus_svv;
    T ut_h = dj * rut + dj_svv * rut_svv;

    shur2[ij] = drdx_local * ur_h +
                drdy_local * us_h +
                drdz_local * ut_h;
    shus2[ij_p] = dsdx_local * ur_h +
                  dsdy_local * us_h +
                  dsdz_local * ut_h;
    rut2 =      dtdx_local * ur_h +
                dtdy_local * us_h +
                dtdz_local * ut_h;

    __syncthreads();

    T uwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur2[l+j*LX] * shdx[l+i*(LX+1)];
      ruw[l] += rut2 * shdz[k+l*(LX+1)];
      uwijke += shus2[i+l*(LX+1)] * shdy[l + j*(LX+1)];
    }
    __syncthreads();
    ruw[k] += uwijke; 
  }
#pragma unroll
  for (int k = 0; k < LX; ++k){
    au[ij + k*LX*LX + ele] = ruw[k];
  }
}

#endif // __MATH_AX_HELM_SVV_KERNEL_H__
