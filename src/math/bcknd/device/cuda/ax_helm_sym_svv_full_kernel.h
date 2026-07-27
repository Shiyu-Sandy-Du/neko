#ifndef __MATH_AX_HELM_SYM_SVV_FULL_KERNEL_H__
#define __MATH_AX_HELM_SYM_SVV_FULL_KERNEL_H__
/*
 Copyright (c) 2026, The Neko Authors
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
 * Fused device kernel for the symmetric full-stress SVV operator.
 *
 * It evaluates div(mu (grad(u) + grad(u)^T)) together with the symmetric
 * trial- and test-filtered SVV contribution. Identity filter matrices disable
 * SVV in inactive directions. No global derivative or stress work arrays are
 * required.
 */
template<typename T, const int LX>
__global__ void ax_helm_sym_svv_full_kernel(
    T * __restrict__ au,
    T * __restrict__ av,
    T * __restrict__ aw,
    const T * __restrict__ u,
    const T * __restrict__ v,
    const T * __restrict__ w,
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
    const T * __restrict__ jacinv,
    const T * __restrict__ w3,
    const T * __restrict__ h1_svv,
    const T * __restrict__ filter_r,
    const T * __restrict__ filter_s,
    const T * __restrict__ filter_t) {

  extern __shared__ T shared[];
  T *shfield = shared;
  T *shwork = shared + LX * LX * LX;

  const int e = blockIdx.x;
  const int i = threadIdx.x;
  const int j = threadIdx.y;
  const int ij = i + j * LX;
  const int lx2 = LX * LX;
  const int elem = e * LX * lx2;

#pragma unroll 1
  for (int k = 0; k < LX; ++k) {
    const int index = ij + k * lx2 + elem;
    au[index] = 0.0;
    av[index] = 0.0;
    aw[index] = 0.0;
  }

  // Process each output component and reference-flux direction.
#pragma unroll 1
  for (int output = 0; output < 3; ++output) {
    T *result = output == 0 ? au : (output == 1 ? av : aw);

#pragma unroll 1
    for (int target = 0; target < 3; ++target) {

      // Ordinary full-stress flux.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        shfield[ij + k * lx2] = 0.0;
      }
      __syncthreads();

#pragma unroll 1
      for (int velocity = 0; velocity < 3; ++velocity) {
        const T *field = velocity == 0 ? u : (velocity == 1 ? v : w);

#pragma unroll 1
        for (int source = 0; source < 3; ++source) {
#pragma unroll 1
          for (int k = 0; k < LX; ++k) {
            T derivative = 0.0;
#pragma unroll
            for (int l = 0; l < LX; ++l) {
              if (source == 0) {
                derivative += dx[i + l * LX] *
                              field[l + j * LX + k * lx2 + elem];
              }
              else if (source == 1) {
                derivative += dy[j + l * LX] *
                              field[i + l * LX + k * lx2 + elem];
              }
              else {
                derivative += dz[k + l * LX] *
                              field[ij + l * lx2 + elem];
              }
            }
            shwork[ij + k * lx2] = derivative;
          }
          __syncthreads();

#pragma unroll 1
          for (int k = 0; k < LX; ++k) {
            const int ijk = ij + k * lx2;
            const int index = ijk + elem;
            T target_x, target_y, target_z;
            T source_x, source_y, source_z;

            if (target == 0) {
              target_x = drdx[index];
              target_y = drdy[index];
              target_z = drdz[index];
            }
            else if (target == 1) {
              target_x = dsdx[index];
              target_y = dsdy[index];
              target_z = dsdz[index];
            }
            else {
              target_x = dtdx[index];
              target_y = dtdy[index];
              target_z = dtdz[index];
            }

            if (source == 0) {
              source_x = drdx[index];
              source_y = drdy[index];
              source_z = drdz[index];
            }
            else if (source == 1) {
              source_x = dsdx[index];
              source_y = dsdy[index];
              source_z = dsdz[index];
            }
            else {
              source_x = dtdx[index];
              source_y = dtdy[index];
              source_z = dtdz[index];
            }

            const T target_velocity =
                velocity == 0 ? target_x :
                (velocity == 1 ? target_y : target_z);
            const T source_output =
                output == 0 ? source_x :
                (output == 1 ? source_y : source_z);
            T coefficient = target_velocity * source_output;
            if (velocity == output) {
              coefficient += target_x * source_x +
                             target_y * source_y +
                             target_z * source_z;
            }
            shfield[ijk] += h1[index] * w3[ijk] * jacinv[index] *
                            coefficient * shwork[ijk];
          }
          __syncthreads();
        }
      }

      // Apply the transposed derivative to the ordinary flux.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        T value = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          if (target == 0) {
            value += dx[l + i * LX] *
                     shfield[l + j * LX + k * lx2];
          }
          else if (target == 1) {
            value += dy[l + j * LX] *
                     shfield[i + l * LX + k * lx2];
          }
          else {
            value += dz[l + k * LX] *
                     shfield[ij + l * lx2];
          }
        }
        result[ij + k * lx2 + elem] += value;
      }
      __syncthreads();

      // Symmetric SVV flux, G Q_hat D u.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        shfield[ij + k * lx2] = 0.0;
      }
      __syncthreads();

#pragma unroll 1
      for (int velocity = 0; velocity < 3; ++velocity) {
        const T *field = velocity == 0 ? u : (velocity == 1 ? v : w);

#pragma unroll 1
        for (int source = 0; source < 3; ++source) {
#pragma unroll 1
          for (int k = 0; k < LX; ++k) {
            T derivative = 0.0;
#pragma unroll
            for (int l = 0; l < LX; ++l) {
              if (source == 0) {
                derivative += dx[i + l * LX] *
                              field[l + j * LX + k * lx2 + elem];
              }
              else if (source == 1) {
                derivative += dy[j + l * LX] *
                              field[i + l * LX + k * lx2 + elem];
              }
              else {
                derivative += dz[k + l * LX] *
                              field[ij + l * lx2 + elem];
              }
            }
            shwork[ij + k * lx2] = derivative;
          }
          __syncthreads();

#pragma unroll 1
          for (int k = 0; k < LX; ++k) {
            const int ijk = ij + k * lx2;
            const int index = ijk + elem;
            T filtered = 0.0;
            T target_x, target_y, target_z;
            T source_x, source_y, source_z;

#pragma unroll
            for (int l = 0; l < LX; ++l) {
              if (source == 0) {
                filtered += filter_r[i + l * LX] *
                            shwork[l + j * LX + k * lx2];
              }
              else if (source == 1) {
                filtered += filter_s[l + j * LX] *
                            shwork[i + l * LX + k * lx2];
              }
              else {
                filtered += filter_t[l + k * LX] *
                            shwork[ij + l * lx2];
              }
            }

            if (target == 0) {
              target_x = drdx[index];
              target_y = drdy[index];
              target_z = drdz[index];
            }
            else if (target == 1) {
              target_x = dsdx[index];
              target_y = dsdy[index];
              target_z = dsdz[index];
            }
            else {
              target_x = dtdx[index];
              target_y = dtdy[index];
              target_z = dtdz[index];
            }

            if (source == 0) {
              source_x = drdx[index];
              source_y = drdy[index];
              source_z = drdz[index];
            }
            else if (source == 1) {
              source_x = dsdx[index];
              source_y = dsdy[index];
              source_z = dsdz[index];
            }
            else {
              source_x = dtdx[index];
              source_y = dtdy[index];
              source_z = dtdz[index];
            }

            const T target_velocity =
                velocity == 0 ? target_x :
                (velocity == 1 ? target_y : target_z);
            const T source_output =
                output == 0 ? source_x :
                (output == 1 ? source_y : source_z);
            T coefficient = target_velocity * source_output;
            if (velocity == output) {
              coefficient += target_x * source_x +
                             target_y * source_y +
                             target_z * source_z;
            }
            shfield[ijk] += h1_svv[index] * w3[ijk] * jacinv[index] *
                            coefficient * (shwork[ijk] - filtered);
          }
          __syncthreads();
        }
      }

      // Test-function-side complementary filter.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        const int ijk = ij + k * lx2;
        T filtered = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          if (target == 0) {
            filtered += filter_r[l + i * LX] *
                        shfield[l + j * LX + k * lx2];
          }
          else if (target == 1) {
            filtered += filter_s[j + l * LX] *
                        shfield[i + l * LX + k * lx2];
          }
          else {
            filtered += filter_t[k + l * LX] *
                        shfield[ij + l * lx2];
          }
        }
        shwork[ijk] = shfield[ijk] - filtered;
      }
      __syncthreads();

      // Apply the transposed derivative to the symmetric SVV flux.
#pragma unroll 1
      for (int k = 0; k < LX; ++k) {
        T value = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          if (target == 0) {
            value += dx[l + i * LX] *
                     shwork[l + j * LX + k * lx2];
          }
          else if (target == 1) {
            value += dy[l + j * LX] *
                     shwork[i + l * LX + k * lx2];
          }
          else {
            value += dz[l + k * LX] *
                     shwork[ij + l * lx2];
          }
        }
        result[ij + k * lx2 + elem] += value;
      }
      __syncthreads();
    }
  }
}

#endif
