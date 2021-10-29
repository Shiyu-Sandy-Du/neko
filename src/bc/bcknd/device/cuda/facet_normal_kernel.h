/**
 * Computes the linear index for area and normal arrays
 * @note Fortran indexing input, C indexing output
 */

#define coef_normal_area_idx(i, j, k, l, lx, nf) \
  (((i) + (lx) * (((j) - 1) + (lx) * (((k) - 1) + (nf) * (((l) - 1))))) - 1)

/**
 * Device function to compute i,j,k,e indices from a linear index
 * @note Assumes idx is a Fortran index 
 */
__device__
void nonlinear_index(const int idx, const int lx, int *index) {
  index[3] = idx/(lx * lx * lx);
  index[2] = (idx - (lx*lx*lx)*index[3])/(lx * lx);
  index[1] = (idx - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) / lx;
  index[0] = (idx - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) - lx*index[1];
  index[1]++;
  index[2]++;
  index[3]++;
}


/**
 * Device kernel for vector apply for a symmetry condition
 */
template< typename T >
__global__
void facet_normal_apply_surfvec_kernel(const int * __restrict__ msk,
				       const int * __restrict__ facet,
				       T * __restrict__ x,
				       T * __restrict__ y,
				       T * __restrict__ z,
				       const T * __restrict__ u,
				       const T * __restrict__ v,
				       const T * __restrict__ w,
				       const T * __restrict__ nx,
				       const T * __restrict__ ny,
				       const T * __restrict__ nz,
				       const T * __restrict__ area,
				       const int lx,
				       const int m) {
  int index[4];
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    const int f = (facet[i]);
    nonlinear_index(msk[i], lx, index);


    switch(f) {
    case 1:
    case 2:
      {
	const int na_idx = coef_normal_area_idx(index[1], index[2],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }
    case 3:
    case 4:
      {
	const int na_idx = coef_normal_area_idx(index[0], index[2],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }
    case 5:
    case 6:
      {
	const int na_idx = coef_normal_area_idx(index[0], index[1],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }    
    }
  }
}
