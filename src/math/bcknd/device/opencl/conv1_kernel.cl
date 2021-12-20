/**
 * Device kernel for convective terms
 */

#define DEFINE_CONV1_KERNEL(LX, CHUNKS)                                        \
__kernel void conv1_kernel_lx##LX(__global real * __restrict__ du,             \
				  __global const real * __restrict__ u,        \
				  __global const real * __restrict__ vx,       \
				  __global const real * __restrict__ vy,       \
				  __global const real * __restrict__ vz,       \
				  __global const real * __restrict__ dx,       \
				  __global const real * __restrict__ dy,       \
				  __global const real * __restrict__ dz,       \
				  __global const real * __restrict__ drdx,     \
				  __global const real * __restrict__ dsdx,     \
				  __global const real * __restrict__ dtdx,     \
				  __global const real * __restrict__ drdy,     \
				  __global const real * __restrict__ dsdy,     \
				  __global const real * __restrict__ dtdy,     \
				  __global const real * __restrict__ drdz,     \
				  __global const real * __restrict__ dsdz,     \
				  __global const real * __restrict__ dtdz,     \
				  __global const real * __restrict__ jacinv) { \
                                                                               \
  __local real shu[LX * LX * LX];                                              \
                                                                               \
  __local real shvx[LX * LX * LX];                                             \
  __local real shvy[LX * LX * LX];                                             \
  __local real shvz[LX * LX * LX];                                             \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shjacinv[LX * LX * LX];                                         \
                                                                               \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[j] = u[j + e * LX * LX * LX];                                          \
                                                                               \
    shvx[j] = vx[j + e * LX * LX * LX];                                        \
    shvy[j] = vy[j + e * LX * LX * LX];                                        \
    shvz[j] = vz[j + e * LX * LX * LX];                                        \
                                                                               \
    shjacinv[j] = jacinv[j + e * LX * LX * LX];                                \
                                                                               \
    j = j + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int n = 0; n < nchunks; n++) {                                          \
    const int ijk = iii + n * CHUNKS;                                          \
    const int jk = ijk / LX;                                                   \
    i = ijk - jk * LX;                                                         \
    k = jk / LX;                                                               \
    j = jk - k * LX;                                                           \
    if ( i < LX && j < LX && k < LX) {                                         \
      real rtmp = 0.0;                                                         \
      real stmp = 0.0;                                                         \
      real ttmp = 0.0;                                                         \
      for (int l = 0; l < LX; l++) {		                               \
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	       \
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];              \
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];              \
      }                                                                        \
                                                                               \
      du[ijk + e * LX * LX * LX] = shjacinv[ijk] *                             \
	(shvx[ijk] * (drdx[ijk + e * LX * LX * LX] * rtmp                      \
		      + dsdx[ijk + e * LX * LX * LX] * stmp                    \
		      + dtdx[ijk + e * LX * LX * LX] * ttmp)                   \
	 + shvy[ijk] * (drdy[ijk + e * LX * LX * LX] * rtmp                    \
			+ dsdy[ijk + e * LX * LX * LX] * stmp                  \
			+ dtdy[ijk + e * LX * LX * LX] * ttmp)                 \
	 + shvz[ijk] * (drdz[ijk + e * LX * LX * LX] * rtmp                    \
			+ dsdz[ijk + e * LX * LX * LX] * stmp                  \
			+ dtdz[ijk + e * LX * LX * LX] * ttmp));               \
    }                                                                          \
  }                                                                            \
}                                                                              

DEFINE_CONV1_KERNEL(2, 256)
DEFINE_CONV1_KERNEL(3, 256)
DEFINE_CONV1_KERNEL(4, 256)
DEFINE_CONV1_KERNEL(5, 256)
DEFINE_CONV1_KERNEL(6, 256)
DEFINE_CONV1_KERNEL(7, 256)
DEFINE_CONV1_KERNEL(8, 256)
DEFINE_CONV1_KERNEL(9, 256)

