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
#include <device/device_config.h>
#include <device/cuda/check.h>
#include "ax_helm_svv_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {

  /**
   * Fortran wrapper for device CUDA Ax_svv version, part 1 (directional)
   */
  void cuda_ax_helm_svv_part1_rst(void *ur, void *us, void *ut,
                        void *u,
                        void *dx, void *dy, void *dz,
                        void *drdx, void *drdy, void *drdz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_rst_KSTEP(LX)                                                       \
    ax_helm_svv_part1_rst_kernel_kstep<real, LX>                                        \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dy, (real *) dz,              \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rst_KSTEP_PADDED(LX)                                                \
    ax_helm_svv_part1_rst_kernel_kstep_padded<real, LX>                                \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dy, (real *) dz,              \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rst(LX)                                                             \
    case LX:                                                                           \
      CASE_part1_rst_KSTEP(LX);                                                        \
       break

#define CASE_part1_rst_PADDED(LX)                                                      \
    case LX:                                                                           \
      CASE_part1_rst_KSTEP_PADDED(LX);                                                 \
       break

    switch(*lx) {
      CASE_part1_rst(2);
      CASE_part1_rst(3);
      CASE_part1_rst_PADDED(4);
      CASE_part1_rst(5);
      CASE_part1_rst(6);
      CASE_part1_rst(7);
      CASE_part1_rst_PADDED(8);
      CASE_part1_rst(9);
      CASE_part1_rst(10);
      CASE_part1_rst(11);
      CASE_part1_rst(12);
      CASE_part1_rst(13);
      CASE_part1_rst(14);
      CASE_part1_rst(15);
      CASE_part1_rst_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_rs(void *ur, void *us, void *ut,
                        void *u,
                        void *dx, void *dy,
                        void *drdx, void *drdy, void *drdz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_rs_KSTEP(LX)                                                        \
    ax_helm_svv_part1_rs_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dy,                           \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rs_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part1_rs_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dy,                           \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rs(LX)                                                              \
    case LX:                                                                           \
      CASE_part1_rs_KSTEP(LX);                                                         \
       break

#define CASE_part1_rs_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part1_rs_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part1_rs(2);
      CASE_part1_rs(3);
      CASE_part1_rs_PADDED(4);
      CASE_part1_rs(5);
      CASE_part1_rs(6);
      CASE_part1_rs(7);
      CASE_part1_rs_PADDED(8);
      CASE_part1_rs(9);
      CASE_part1_rs(10);
      CASE_part1_rs(11);
      CASE_part1_rs(12);
      CASE_part1_rs(13);
      CASE_part1_rs(14);
      CASE_part1_rs(15);
      CASE_part1_rs_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_rt(void *ur, void *us, void *ut,
                        void *u,
                        void *dx, void *dz,
                        void *drdx, void *drdy, void *drdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_rt_KSTEP(LX)                                                        \
    ax_helm_svv_part1_rt_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dz,                           \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rt_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part1_rt_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx, (real *) dz,                           \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_rt(LX)                                                              \
    case LX:                                                                           \
      CASE_part1_rt_KSTEP(LX);                                                         \
       break

#define CASE_part1_rt_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part1_rt_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part1_rt(2);
      CASE_part1_rt(3);
      CASE_part1_rt_PADDED(4);
      CASE_part1_rt(5);
      CASE_part1_rt(6);
      CASE_part1_rt(7);
      CASE_part1_rt_PADDED(8);
      CASE_part1_rt(9);
      CASE_part1_rt(10);
      CASE_part1_rt(11);
      CASE_part1_rt(12);
      CASE_part1_rt(13);
      CASE_part1_rt(14);
      CASE_part1_rt(15);
      CASE_part1_rt_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_st(void *ur, void *us, void *ut,
                        void *u,
                        void *dy, void *dz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_st_KSTEP(LX)                                                        \
    ax_helm_svv_part1_st_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dy, (real *) dz,                           \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_st_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part1_st_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dy, (real *) dz,                           \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_st(LX)                                                              \
    case LX:                                                                           \
      CASE_part1_st_KSTEP(LX);                                                         \
       break

#define CASE_part1_st_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part1_st_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part1_st(2);
      CASE_part1_st(3);
      CASE_part1_st_PADDED(4);
      CASE_part1_st(5);
      CASE_part1_st(6);
      CASE_part1_st(7);
      CASE_part1_st_PADDED(8);
      CASE_part1_st(9);
      CASE_part1_st(10);
      CASE_part1_st(11);
      CASE_part1_st(12);
      CASE_part1_st(13);
      CASE_part1_st(14);
      CASE_part1_st(15);
      CASE_part1_st_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_r(void *ur, void *us, void *ut,
                        void *u,
                        void *dx,
                        void *drdx, void *drdy, void *drdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_r_KSTEP(LX)                                                         \
    ax_helm_svv_part1_r_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_r_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part1_r_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dx,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_r(LX)                                                               \
    case LX:                                                                           \
      CASE_part1_r_KSTEP(LX);                                                          \
       break

#define CASE_part1_r_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part1_r_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part1_r(2);
      CASE_part1_r(3);
      CASE_part1_r_PADDED(4);
      CASE_part1_r(5);
      CASE_part1_r(6);
      CASE_part1_r(7);
      CASE_part1_r_PADDED(8);
      CASE_part1_r(9);
      CASE_part1_r(10);
      CASE_part1_r(11);
      CASE_part1_r(12);
      CASE_part1_r(13);
      CASE_part1_r(14);
      CASE_part1_r(15);
      CASE_part1_r_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_s(void *ur, void *us, void *ut,
                        void *u,
                        void *dy,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_s_KSTEP(LX)                                                         \
    ax_helm_svv_part1_s_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dy,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_s_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part1_s_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dy,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_s(LX)                                                               \
    case LX:                                                                           \
      CASE_part1_s_KSTEP(LX);                                                          \
       break

#define CASE_part1_s_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part1_s_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part1_s(2);
      CASE_part1_s(3);
      CASE_part1_s_PADDED(4);
      CASE_part1_s(5);
      CASE_part1_s(6);
      CASE_part1_s(7);
      CASE_part1_s_PADDED(8);
      CASE_part1_s(9);
      CASE_part1_s(10);
      CASE_part1_s(11);
      CASE_part1_s(12);
      CASE_part1_s(13);
      CASE_part1_s(14);
      CASE_part1_s(15);
      CASE_part1_s_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part1_t(void *ur, void *us, void *ut,
                        void *u,
                        void *dz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *jacinv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part1_t_KSTEP(LX)                                                         \
    ax_helm_svv_part1_t_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dz,                                        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_t_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part1_t_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) ur, (real *) us, (real *) ut,              \
                                     (real *) u,                                         \
                                     (real *) dz,                                        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) jacinv);                                   \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part1_t(LX)                                                               \
    case LX:                                                                           \
      CASE_part1_t_KSTEP(LX);                                                          \
       break

#define CASE_part1_t_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part1_t_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part1_t(2);
      CASE_part1_t(3);
      CASE_part1_t_PADDED(4);
      CASE_part1_t(5);
      CASE_part1_t(6);
      CASE_part1_t(7);
      CASE_part1_t_PADDED(8);
      CASE_part1_t(9);
      CASE_part1_t(10);
      CASE_part1_t(11);
      CASE_part1_t(12);
      CASE_part1_t(13);
      CASE_part1_t(14);
      CASE_part1_t(15);
      CASE_part1_t_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  /**
   * Fortran wrapper for device CUDA Ax_svv version, part 2 (directional)
   */
  void cuda_ax_helm_svv_part2_rst(void *au,
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
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_rst_KSTEP(LX)                                                       \
    ax_helm_svv_part2_rst_kernel_kstep<real, LX>                                        \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dy, (real *) dz,              \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rst_KSTEP_PADDED(LX)                                                \
    ax_helm_svv_part2_rst_kernel_kstep_padded<real, LX>                                \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dy, (real *) dz,              \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rst(LX)                                                             \
    case LX:                                                                           \
      CASE_part2_rst_KSTEP(LX);                                                        \
       break

#define CASE_part2_rst_PADDED(LX)                                                      \
    case LX:                                                                           \
      CASE_part2_rst_KSTEP_PADDED(LX);                                                 \
       break

    switch(*lx) {
      CASE_part2_rst(2);
      CASE_part2_rst(3);
      CASE_part2_rst_PADDED(4);
      CASE_part2_rst(5);
      CASE_part2_rst(6);
      CASE_part2_rst(7);
      CASE_part2_rst_PADDED(8);
      CASE_part2_rst(9);
      CASE_part2_rst(10);
      CASE_part2_rst(11);
      CASE_part2_rst(12);
      CASE_part2_rst(13);
      CASE_part2_rst(14);
      CASE_part2_rst(15);
      CASE_part2_rst_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_rs(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dx, void *dy,
                        void *h1,
                        void *drdx, void *drdy, void *drdz,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_rs_KSTEP(LX)                                                        \
    ax_helm_svv_part2_rs_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dy,                           \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rs_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part2_rs_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dy,                           \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rs(LX)                                                              \
    case LX:                                                                           \
      CASE_part2_rs_KSTEP(LX);                                                         \
       break

#define CASE_part2_rs_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part2_rs_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part2_rs(2);
      CASE_part2_rs(3);
      CASE_part2_rs_PADDED(4);
      CASE_part2_rs(5);
      CASE_part2_rs(6);
      CASE_part2_rs(7);
      CASE_part2_rs_PADDED(8);
      CASE_part2_rs(9);
      CASE_part2_rs(10);
      CASE_part2_rs(11);
      CASE_part2_rs(12);
      CASE_part2_rs(13);
      CASE_part2_rs(14);
      CASE_part2_rs(15);
      CASE_part2_rs_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_rt(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dx, void *dz,
                        void *h1,
                        void *drdx, void *drdy, void *drdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_rt_KSTEP(LX)                                                        \
    ax_helm_svv_part2_rt_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dz,                           \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rt_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part2_rt_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx, (real *) dz,                           \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_rt(LX)                                                              \
    case LX:                                                                           \
      CASE_part2_rt_KSTEP(LX);                                                         \
       break

#define CASE_part2_rt_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part2_rt_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part2_rt(2);
      CASE_part2_rt(3);
      CASE_part2_rt_PADDED(4);
      CASE_part2_rt(5);
      CASE_part2_rt(6);
      CASE_part2_rt(7);
      CASE_part2_rt_PADDED(8);
      CASE_part2_rt(9);
      CASE_part2_rt(10);
      CASE_part2_rt(11);
      CASE_part2_rt(12);
      CASE_part2_rt(13);
      CASE_part2_rt(14);
      CASE_part2_rt(15);
      CASE_part2_rt_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_st(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dy, void *dz,
                        void *h1,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_st_KSTEP(LX)                                                        \
    ax_helm_svv_part2_st_kernel_kstep<real, LX>                                         \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dy, (real *) dz,                           \
                                     (real *) h1,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_st_KSTEP_PADDED(LX)                                                 \
    ax_helm_svv_part2_st_kernel_kstep_padded<real, LX>                                 \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dy, (real *) dz,                           \
                                     (real *) h1,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_st(LX)                                                              \
    case LX:                                                                           \
      CASE_part2_st_KSTEP(LX);                                                         \
       break

#define CASE_part2_st_PADDED(LX)                                                       \
    case LX:                                                                           \
      CASE_part2_st_KSTEP_PADDED(LX);                                                  \
       break

    switch(*lx) {
      CASE_part2_st(2);
      CASE_part2_st(3);
      CASE_part2_st_PADDED(4);
      CASE_part2_st(5);
      CASE_part2_st(6);
      CASE_part2_st(7);
      CASE_part2_st_PADDED(8);
      CASE_part2_st(9);
      CASE_part2_st(10);
      CASE_part2_st(11);
      CASE_part2_st(12);
      CASE_part2_st(13);
      CASE_part2_st(14);
      CASE_part2_st(15);
      CASE_part2_st_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_r(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dx,
                        void *h1,
                        void *drdx, void *drdy, void *drdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_r_KSTEP(LX)                                                         \
    ax_helm_svv_part2_r_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx,                                        \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_r_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part2_r_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dx,                                        \
                                     (real *) h1,                                        \
                                     (real *) drdx, (real *) drdy, (real *) drdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_r(LX)                                                               \
    case LX:                                                                           \
      CASE_part2_r_KSTEP(LX);                                                          \
       break

#define CASE_part2_r_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part2_r_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part2_r(2);
      CASE_part2_r(3);
      CASE_part2_r_PADDED(4);
      CASE_part2_r(5);
      CASE_part2_r(6);
      CASE_part2_r(7);
      CASE_part2_r_PADDED(8);
      CASE_part2_r(9);
      CASE_part2_r(10);
      CASE_part2_r(11);
      CASE_part2_r(12);
      CASE_part2_r(13);
      CASE_part2_r(14);
      CASE_part2_r(15);
      CASE_part2_r_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_s(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dy,
                        void *h1,
                        void *dsdx, void *dsdy, void *dsdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_s_KSTEP(LX)                                                         \
    ax_helm_svv_part2_s_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dy,                                        \
                                     (real *) h1,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_s_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part2_s_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dy,                                        \
                                     (real *) h1,                                        \
                                     (real *) dsdx, (real *) dsdy, (real *) dsdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_s(LX)                                                               \
    case LX:                                                                           \
      CASE_part2_s_KSTEP(LX);                                                          \
       break

#define CASE_part2_s_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part2_s_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part2_s(2);
      CASE_part2_s(3);
      CASE_part2_s_PADDED(4);
      CASE_part2_s(5);
      CASE_part2_s(6);
      CASE_part2_s(7);
      CASE_part2_s_PADDED(8);
      CASE_part2_s(9);
      CASE_part2_s(10);
      CASE_part2_s(11);
      CASE_part2_s(12);
      CASE_part2_s(13);
      CASE_part2_s(14);
      CASE_part2_s(15);
      CASE_part2_s_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  void cuda_ax_helm_svv_part2_t(void *au,
                        void *ur, void *us, void *ut,
                        void *ur_svv, void *us_svv, void *ut_svv,
                        void *dz,
                        void *h1,
                        void *dtdx, void *dtdy, void *dtdz,
                        void *w3, void *h1_svv, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_part2_t_KSTEP(LX)                                                         \
    ax_helm_svv_part2_t_kernel_kstep<real, LX>                                          \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dz,                                        \
                                     (real *) h1,                                        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_t_KSTEP_PADDED(LX)                                                  \
    ax_helm_svv_part2_t_kernel_kstep_padded<real, LX>                                  \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au,                                        \
                                     (real *) ur, (real *) us, (real *) ut,              \
                                     (real *) ur_svv, (real *) us_svv,                   \
                                     (real *) ut_svv,                                    \
                                     (real *) dz,                                        \
                                     (real *) h1,                                        \
                                     (real *) dtdx, (real *) dtdy, (real *) dtdz,        \
                                     (real *) w3, (real *) h1_svv);                      \
    CUDA_CHECK(cudaGetLastError());

#define CASE_part2_t(LX)                                                               \
    case LX:                                                                           \
      CASE_part2_t_KSTEP(LX);                                                          \
       break

#define CASE_part2_t_PADDED(LX)                                                        \
    case LX:                                                                           \
      CASE_part2_t_KSTEP_PADDED(LX);                                                   \
       break

    switch(*lx) {
      CASE_part2_t(2);
      CASE_part2_t(3);
      CASE_part2_t_PADDED(4);
      CASE_part2_t(5);
      CASE_part2_t(6);
      CASE_part2_t(7);
      CASE_part2_t_PADDED(8);
      CASE_part2_t(9);
      CASE_part2_t(10);
      CASE_part2_t(11);
      CASE_part2_t(12);
      CASE_part2_t(13);
      CASE_part2_t(14);
      CASE_part2_t(15);
      CASE_part2_t_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

}
