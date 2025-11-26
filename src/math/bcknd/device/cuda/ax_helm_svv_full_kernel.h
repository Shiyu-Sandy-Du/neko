#ifndef __MATH_AX_HELM_SVV_FULL_KERNEL_H__
#define __MATH_AX_HELM_SVV_FULL_KERNEL_H__
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
 * Device kernels for Ax helm svv full
 */

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_full_part1_kernel_vector_kstep(T * __restrict__ s11,
                                      T * __restrict__ s22,
                                      T * __restrict__ s33,
                                      T * __restrict__ s12,
                                      T * __restrict__ s13,
                                      T * __restrict__ s23,
                                      const T * __restrict__ u,
                                      const T * __restrict__ v,
                                      const T * __restrict__ w,                            
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
  __shared__ T shv[LX * LX];
  __shared__ T shw[LX * LX];

  T ru[LX];
  T rv[LX];
  T rw[LX];

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
    rv[k] = v[ij + k*LX*LX + ele];
    rw[k] = w[ij + k*LX*LX + ele];
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
    T vttmp = 0.0;
    T wttmp = 0.0;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];

    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*LX] * ru[l];
      vttmp += shdz[k+l*LX] * rv[l];
      wttmp += shdz[k+l*LX] * rw[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

    T vrtmp = 0.0;
    T vstmp = 0.0;

    T wrtmp = 0.0;
    T wstmp = 0.0;

#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*LX] * shu[l+j*LX];
      ustmp += shdy[j+l*LX] * shu[i+l*LX];

      vrtmp += shdx[i+l*LX] * shv[l+j*LX];
      vstmp += shdy[j+l*LX] * shv[i+l*LX];

      wrtmp += shdx[i+l*LX] * shw[l+j*LX];
      wstmp += shdy[j+l*LX] * shw[i+l*LX];
    }
    
    T u1 = 0.0;
    T u2 = 0.0;
    T u3 = 0.0;
    T v1 = 0.0;
    T v2 = 0.0;
    T v3 = 0.0;
    T w1 = 0.0;
    T w2 = 0.0;
    T w3 = 0.0;

    u1 = urtmp * drdx_local + 
         ustmp * dsdx_local + 
         uttmp * dtdx_local;
    u2 = urtmp * drdy_local + 
         ustmp * dsdy_local + 
         uttmp * dtdy_local;
    u3 = urtmp * drdz_local + 
         ustmp * dsdz_local + 
         uttmp * dtdz_local;

    v1 = vrtmp * drdx_local + 
         vstmp * dsdx_local + 
         vttmp * dtdx_local;
    v2 = vrtmp * drdy_local + 
         vstmp * dsdy_local + 
         vttmp * dtdy_local;
    v3 = vrtmp * drdz_local + 
         vstmp * dsdz_local + 
         vttmp * dtdz_local;

    w1 = wrtmp * drdx_local + 
         wstmp * dsdx_local + 
         wttmp * dtdx_local;
    w2 = wrtmp * drdy_local + 
         wstmp * dsdy_local + 
         wttmp * dtdy_local;
    w3 = wrtmp * drdz_local + 
         wstmp * dsdz_local + 
         wttmp * dtdz_local;

    s11[ij + k*LX*LX + ele] = dj*(u1 + u1);
    s22[ij + k*LX*LX + ele] = dj*(v2 + v2);
    s33[ij + k*LX*LX + ele] = dj*(w3 + w3);
    s12[ij + k*LX*LX + ele] = dj*(u2 + v1);
    s13[ij + k*LX*LX + ele] = dj*(u3 + w1);
    s23[ij + k*LX*LX + ele] = dj*(v3 + w2);   
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_full_part1_kernel_vector_kstep_padded(T * __restrict__ s11,
                                             T * __restrict__ s22,
                                             T * __restrict__ s33,
                                             T * __restrict__ s12,
                                             T * __restrict__ s13,
                                             T * __restrict__ s23,
                                             const T * __restrict__ u,
                                             const T * __restrict__ v,
                                             const T * __restrict__ w,                            
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
  __shared__ T shv[LX * (LX+1)];
  __shared__ T shw[LX * (LX+1)];

  T ru[LX];
  T rv[LX];
  T rw[LX];

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
    rv[k] = v[ij + k*LX*LX + ele];
    rw[k] = w[ij + k*LX*LX + ele];
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
    T vttmp = 0.0;
    T wttmp = 0.0;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];

    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*(LX+1)] * ru[l];
      vttmp += shdz[k+l*(LX+1)] * rv[l];
      wttmp += shdz[k+l*(LX+1)] * rw[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

    T vrtmp = 0.0;
    T vstmp = 0.0;

    T wrtmp = 0.0;
    T wstmp = 0.0;

#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*(LX+1)] * shu[l+j*(LX+1)];
      ustmp += shdy[j+l*(LX+1)] * shu[i+l*(LX+1)];

      vrtmp += shdx[i+l*(LX+1)] * shv[l+j*(LX+1)];
      vstmp += shdy[j+l*(LX+1)] * shv[i+l*(LX+1)];

      wrtmp += shdx[i+l*(LX+1)] * shw[l+j*(LX+1)];
      wstmp += shdy[j+l*(LX+1)] * shw[i+l*(LX+1)];
    }

    T u1 = 0.0;
    T u2 = 0.0;
    T u3 = 0.0;
    T v1 = 0.0;
    T v2 = 0.0;
    T v3 = 0.0;
    T w1 = 0.0;
    T w2 = 0.0;
    T w3 = 0.0;

    u1 = urtmp * drdx_local + 
         ustmp * dsdx_local + 
         uttmp * dtdx_local;
    u2 = urtmp * drdy_local + 
         ustmp * dsdy_local + 
         uttmp * dtdy_local;
    u3 = urtmp * drdz_local + 
         ustmp * dsdz_local + 
         uttmp * dtdz_local;

    v1 = vrtmp * drdx_local + 
         vstmp * dsdx_local + 
         vttmp * dtdx_local;
    v2 = vrtmp * drdy_local + 
         vstmp * dsdy_local + 
         vttmp * dtdy_local;
    v3 = vrtmp * drdz_local + 
         vstmp * dsdz_local + 
         vttmp * dtdz_local;

    w1 = wrtmp * drdx_local + 
         wstmp * dsdx_local + 
         wttmp * dtdx_local;
    w2 = wrtmp * drdy_local + 
         wstmp * dsdy_local + 
         wttmp * dtdy_local;
    w3 = wrtmp * drdz_local + 
         wstmp * dsdz_local + 
         wttmp * dtdz_local;

    s11[ij + k*LX*LX + ele] = dj*(u1 + u1);
    s22[ij + k*LX*LX + ele] = dj*(v2 + v2);
    s33[ij + k*LX*LX + ele] = dj*(w3 + w3);
    s12[ij + k*LX*LX + ele] = dj*(u2 + v1);
    s13[ij + k*LX*LX + ele] = dj*(u3 + w1);
    s23[ij + k*LX*LX + ele] = dj*(v3 + w2);
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_full_part2_kernel_vector_kstep(T * __restrict__ au,
                                 T * __restrict__ av,
                                 T * __restrict__ aw,
                                 const T * __restrict__ s11,
                                 const T * __restrict__ s22,
                                 const T * __restrict__ s33,
                                 const T * __restrict__ s12,
                                 const T * __restrict__ s13,
                                 const T * __restrict__ s23,
                                 const T * __restrict__ s11_svv,
                                 const T * __restrict__ s22_svv,
                                 const T * __restrict__ s33_svv,
                                 const T * __restrict__ s12_svv,
                                 const T * __restrict__ s13_svv,
                                 const T * __restrict__ s23_svv,
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
  __shared__ T shvr2[LX * LX];
  __shared__ T shvs2[LX * LX];
  T rvt2;
  __shared__ T shwr2[LX * LX];
  __shared__ T shws2[LX * LX];
  T rwt2;
  
  T rs11[LX];
  T rs22[LX];
  T rs33[LX];
  T rs12[LX];
  T rs13[LX];
  T rs23[LX];  

  T rs11_svv[LX];
  T rs22_svv[LX];
  T rs33_svv[LX];
  T rs12_svv[LX];
  T rs13_svv[LX];
  T rs23_svv[LX];

  T ruw[LX];
  T rvw[LX];
  T rww[LX];

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
    rs11[k] = s11[ij + k*LX*LX + ele];
    rs22[k] = s22[ij + k*LX*LX + ele];
    rs33[k] = s33[ij + k*LX*LX + ele];
    rs12[k] = s12[ij + k*LX*LX + ele];
    rs13[k] = s13[ij + k*LX*LX + ele];
    rs23[k] = s23[ij + k*LX*LX + ele];
    rs11_svv[k] = s11[ij + k*LX*LX + ele] - s11_svv[ij + k*LX*LX + ele];
    rs22_svv[k] = s22[ij + k*LX*LX + ele] - s22_svv[ij + k*LX*LX + ele];
    rs33_svv[k] = s33[ij + k*LX*LX + ele] - s33_svv[ij + k*LX*LX + ele];
    rs12_svv[k] = s12[ij + k*LX*LX + ele] - s12_svv[ij + k*LX*LX + ele];
    rs13_svv[k] = s13[ij + k*LX*LX + ele] - s13_svv[ij + k*LX*LX + ele];
    rs23_svv[k] = s23[ij + k*LX*LX + ele] - s23_svv[ij + k*LX*LX + ele];
    ruw[k] = 0.0;
    rvw[k] = 0.0;
    rww[k] = 0.0;
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

    T rs11_h = dj * rs11[k] + dj_svv * rs11_svv[k];
    T rs22_h = dj * rs22[k] + dj_svv * rs22_svv[k];
    T rs33_h = dj * rs33[k] + dj_svv * rs33_svv[k];
    T rs12_h = dj * rs12[k] + dj_svv * rs12_svv[k];
    T rs13_h = dj * rs13[k] + dj_svv * rs13_svv[k];
    T rs23_h = dj * rs23[k] + dj_svv * rs23_svv[k];

    shur2[ij] = drdx_local * rs11_h +
                drdy_local * rs12_h +
                drdz_local * rs13_h;
    shus2[ij] = dsdx_local * rs11_h +
                dsdy_local * rs12_h +
                dsdz_local * rs13_h;
    rut2 =      dtdx_local * rs11_h +
                dtdy_local * rs12_h +
                dtdz_local * rs13_h;
    shvr2[ij] = drdx_local * rs12_h +
                drdy_local * rs22_h +
                drdz_local * rs23_h;
    shvs2[ij] = dsdx_local * rs12_h +
                dsdy_local * rs22_h +
                dsdz_local * rs23_h;
    rvt2 =      dtdx_local * rs12_h +
                dtdy_local * rs22_h +
                dtdz_local * rs23_h;
    shwr2[ij] = drdx_local * rs13_h +
                drdy_local * rs23_h +
                drdz_local * rs33_h;
    shws2[ij] = dsdx_local * rs13_h +
                dsdy_local * rs23_h +
                dsdz_local * rs33_h;
    rwt2 =      dtdx_local * rs13_h +
                dtdy_local * rs23_h +
                dtdz_local * rs33_h;

    __syncthreads();

    T uwijke = 0.0;
    T vwijke = 0.0;
    T wwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur2[l+j*LX] * shdx[l+i*LX];
      ruw[l] += rut2 * shdz[k+l*LX];
      uwijke += shus2[i+l*LX] * shdy[l + j*LX];
      vwijke += shvr2[l+j*LX] * shdx[l+i*LX];
      rvw[l] += rvt2 * shdz[k+l*LX];
      vwijke += shvs2[i+l*LX] * shdy[l + j*LX];
      wwijke += shwr2[l+j*LX] * shdx[l+i*LX];
      rww[l] += rwt2 * shdz[k+l*LX];
      wwijke += shws2[i+l*LX] * shdy[l + j*LX];
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }
#pragma unroll
  for (int k = 0; k < LX; ++k){
    au[ij + k*LX*LX + ele] = ruw[k];
    av[ij + k*LX*LX + ele] = rvw[k];
    aw[ij + k*LX*LX + ele] = rww[k];
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  ax_helm_svv_full_part2_kernel_vector_kstep_padded(T * __restrict__ au,
                                 T * __restrict__ av,
                                 T * __restrict__ aw,
                                 const T * __restrict__ s11,
                                 const T * __restrict__ s22,
                                 const T * __restrict__ s33,
                                 const T * __restrict__ s12,
                                 const T * __restrict__ s13,
                                 const T * __restrict__ s23,
                                 const T * __restrict__ s11_svv,
                                 const T * __restrict__ s22_svv,
                                 const T * __restrict__ s33_svv,
                                 const T * __restrict__ s12_svv,
                                 const T * __restrict__ s13_svv,
                                 const T * __restrict__ s23_svv,
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
  __shared__ T shvr2[LX * LX];
  __shared__ T shvs2[LX * (LX+1)];
  T rvt2;
  __shared__ T shwr2[LX * LX];
  __shared__ T shws2[LX * (LX+1)];
  T rwt2;
  
  T rs11[LX];
  T rs22[LX];
  T rs33[LX];
  T rs12[LX];
  T rs13[LX];
  T rs23[LX];  

  T rs11_svv[LX];
  T rs22_svv[LX];
  T rs33_svv[LX];
  T rs12_svv[LX];
  T rs13_svv[LX];
  T rs23_svv[LX];

  T ruw[LX];
  T rvw[LX];
  T rww[LX];

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
    rs11[k] = s11[ij + k*LX*LX + ele];
    rs22[k] = s22[ij + k*LX*LX + ele];
    rs33[k] = s33[ij + k*LX*LX + ele];
    rs12[k] = s12[ij + k*LX*LX + ele];
    rs13[k] = s13[ij + k*LX*LX + ele];
    rs23[k] = s23[ij + k*LX*LX + ele];
    rs11_svv[k] = s11[ij + k*LX*LX + ele] - s11_svv[ij + k*LX*LX + ele];
    rs22_svv[k] = s22[ij + k*LX*LX + ele] - s22_svv[ij + k*LX*LX + ele];
    rs33_svv[k] = s33[ij + k*LX*LX + ele] - s33_svv[ij + k*LX*LX + ele];
    rs12_svv[k] = s12[ij + k*LX*LX + ele] - s12_svv[ij + k*LX*LX + ele];
    rs13_svv[k] = s13[ij + k*LX*LX + ele] - s13_svv[ij + k*LX*LX + ele];
    rs23_svv[k] = s23[ij + k*LX*LX + ele] - s23_svv[ij + k*LX*LX + ele];
    ruw[k] = 0.0;
    rvw[k] = 0.0;
    rww[k] = 0.0;
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

    T rs11_h = dj * rs11[k] + dj_svv * rs11_svv[k];
    T rs22_h = dj * rs22[k] + dj_svv * rs22_svv[k];
    T rs33_h = dj * rs33[k] + dj_svv * rs33_svv[k];
    T rs12_h = dj * rs12[k] + dj_svv * rs12_svv[k];
    T rs13_h = dj * rs13[k] + dj_svv * rs13_svv[k];
    T rs23_h = dj * rs23[k] + dj_svv * rs23_svv[k];

    shur2[ij] = drdx_local * rs11_h +
                drdy_local * rs12_h +
                drdz_local * rs13_h;
    shus2[ij_p] = dsdx_local * rs11_h +
                dsdy_local * rs12_h +
                dsdz_local * rs13_h;
    rut2 =      dtdx_local * rs11_h +
                dtdy_local * rs12_h +
                dtdz_local * rs13_h;
    shvr2[ij] = drdx_local * rs12_h +
                drdy_local * rs22_h +
                drdz_local * rs23_h;
    shvs2[ij_p] = dsdx_local * rs12_h +
                dsdy_local * rs22_h +
                dsdz_local * rs23_h;
    rvt2 =      dtdx_local * rs12_h +
                dtdy_local * rs22_h +
                dtdz_local * rs23_h;
    shwr2[ij] = drdx_local * rs13_h +
                drdy_local * rs23_h +
                drdz_local * rs33_h;
    shws2[ij_p] = dsdx_local * rs13_h +
                dsdy_local * rs23_h +
                dsdz_local * rs33_h;
    rwt2 =      dtdx_local * rs13_h +
                dtdy_local * rs23_h +
                dtdz_local * rs33_h;

    __syncthreads();

    T uwijke = 0.0;
    T vwijke = 0.0;
    T wwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur2[l+j*LX] * shdx[l+i*(LX+1)];
      ruw[l] += rut2 * shdz[k+l*(LX+1)];
      uwijke += shus2[i+l*(LX+1)] * shdy[l + j*(LX+1)];

      vwijke += shvr2[l+j*LX] * shdx[l+i*(LX+1)];
      rvw[l] += rvt2 * shdz[k+l*(LX+1)];
      vwijke += shvs2[i+l*(LX+1)] * shdy[l + j*(LX+1)];

      wwijke += shwr2[l+j*LX] * shdx[l+i*(LX+1)];
      rww[l] += rwt2 * shdz[k+l*(LX+1)];
      wwijke += shws2[i+l*(LX+1)] * shdy[l + j*(LX+1)];
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }
#pragma unroll
  for (int k = 0; k < LX; ++k){
    au[ij + k*LX*LX + ele] = ruw[k];
    av[ij + k*LX*LX + ele] = rvw[k];
    aw[ij + k*LX*LX + ele] = rww[k];
  }
}

#endif // __MATH_AX_HELM_SVV_FULL_KERNEL_H__
