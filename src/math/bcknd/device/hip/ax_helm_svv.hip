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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ax_helm_svv_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {

  /**
   * Fortran wrapper for device HIP Ax_svv version, part 1
   */
  void hip_ax_helm_svv_part1(void *ur, void *us, void *ut,
                        void *u,
                        void *dx, void *dy, void *dz,
                        void *drdx, void *drdy, void *drdz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);

#define CASE_part1_KSTEP(LX)                                                        \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(ax_helm_svv_part1_kernel_kstep<real, LX> ), \
                        nblcks, nthrds, 0,                                          \
                        (hipStream_t) glb_cmd_queue,                                \
                        (real *) ur, (real *) us, (real *) ut,                      \
                        (real *) u,                                                 \
                        (real *) dx, (real *) dy, (real *) dz,                      \
                        (real *) drdx, (real *) drdy, (real *) drdz,                \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,                \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,                \
                        (real *) jacinv);                                           \
    HIP_CHECK(hipGetLastError());

#define CASE_part1_KSTEP_PADDED(LX)                                                        \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(ax_helm_svv_part1_kernel_kstep_padded<real, LX> ), \
                        nblcks, nthrds, 0,                                                 \
                        (hipStream_t) glb_cmd_queue,                                       \
                        (real *) ur, (real *) us, (real *) ut,                             \
                        (real *) u,                                                        \
                        (real *) dx, (real *) dy, (real *) dz,                             \
                        (real *) drdx, (real *) drdy, (real *) drdz,                       \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,                       \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,                       \
                        (real *) jacinv);                                                  \
    HIP_CHECK(hipGetLastError());

#define CASE_part1(LX)                                                          \
    case LX:                                                                    \
      CASE_part1_KSTEP(LX);                                                     \
       break

#define CASE_part1_PADDED(LX)                                                   \
    case LX:                                                                    \
      CASE_part1_KSTEP_PADDED(LX);                                              \
       break

    switch(*lx) {
      CASE_part1(2);
      CASE_part1(3);
      CASE_part1_PADDED(4);
      CASE_part1(5);
      CASE_part1(6);
      CASE_part1(7);
      CASE_part1_PADDED(8);
      CASE_part1(9);
      CASE_part1(10);
      CASE_part1(11);
      CASE_part1(12);
      CASE_part1(13);
      CASE_part1(14);
      CASE_part1(15);
      CASE_part1_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  /**
   * Fortran wrapper for device HIP Ax_svv version, part 2
   */
  void hip_ax_helm_svv_part2(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dx, void *dy, void *dz,
                        void *h1,
                        void *drdx, void *drdy, void *drdz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);

#define CASE_part2_KSTEP(LX)                                                         \
    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_svv_part2_kernel_kstep<real, LX> ), \
                        nblcks, nthrds, 0,                                           \
                        (hipStream_t) glb_cmd_queue,                                 \
                        (real *) au,                                                 \
                        (real *) ur, (real *) us, (real *) ut,                       \
                        (real *) ur_svv, (real *) us_svv,                            \
                        (real *) ut_svv,                                             \
                        (real *) dx, (real *) dy, (real *) dz,                       \
                        (real *) h1,                                                 \
                        (real *) drdx, (real *) drdy, (real *) drdz,                 \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,                 \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,                 \
                        (real *) w3, (real *) h1_svv);                               \
    HIP_CHECK(hipGetLastError());


#define CASE_part2_KSTEP_PADDED(LX)                                                         \
    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_svv_part2_kernel_kstep_padded<real, LX> ), \
                        nblcks, nthrds, 0,                                                  \
                        (hipStream_t) glb_cmd_queue,                                        \
                        (real *) au,                                                        \
                        (real *) ur, (real *) us, (real *) ut,                              \
                        (real *) ur_svv, (real *) us_svv,                                   \
                        (real *) ut_svv,                                                    \
                        (real *) dx, (real *) dy, (real *) dz,                              \
                        (real *) h1,                                                        \
                        (real *) drdx, (real *) drdy, (real *) drdz,                        \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,                        \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,                        \
                        (real *) w3, (real *) h1_svv);                                      \
    HIP_CHECK(hipGetLastError());

#define CASE_part2(LX)                                                          \
    case LX:                                                                    \
      CASE_part2_KSTEP(LX);                                                     \
       break

#define CASE_part2_PADDED(LX)                                                   \
    case LX:                                                                    \
      CASE_part2_KSTEP_PADDED(LX);                                              \
       break

    switch(*lx) {
      CASE_part2(2);
      CASE_part2(3);
      CASE_part2_PADDED(4);
      CASE_part2(5);
      CASE_part2(6);
      CASE_part2(7);
      CASE_part2_PADDED(8);
      CASE_part2(9);
      CASE_part2(10);
      CASE_part2(11);
      CASE_part2(12);
      CASE_part2(13);
      CASE_part2(14);
      CASE_part2(15);
      CASE_part2_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

}
