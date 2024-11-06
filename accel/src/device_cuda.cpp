/* -*- c++ -*- */

#if defined(_GPU_CUDA)

#include "device.h"

#include <stdio.h>

// ALGO1: all threads calculate all DF
// ALGO2: threads calculate subset of DF
#define _ALGO2

#if defined(_CUDA_NVTX)
#include "nvToolsExt.h"
#endif

/* ---------------------------------------------------------------------- */

inline __host__ __device__ double pow_int(double a, int n)
{
  double result;
  if(n == 0) result = 1.0;
  else if(n == 1) result = a;
  else if(n == 2) result = a*a;
  else if(n == 3) result = a*a*a;
  else if(n == 4) {
    double a2 = a*a;
    result = a2*a2;
  } else if(n == 5) {
    double a2 = a*a;
    result = a2*a2*a;
  } else if(n == 6) {
    double a2 = a*a;
    result = a2*a2*a2;
  } else if(n == 7) {
    double a2 = a*a;
    result = a2*a2*a2*a;
  } else if(n == 8) {
    double a2 = a*a;
    result = a2*a2*a2*a2;
  } else if(n == 9) {
    double a2 = a*a;
    result = a2*a2*a2*a2*a;
  }

  return result;
}

/* ---------------------------------------------------------------------- */

#ifdef _ALGO2

__global__ void _compute_pip_algo2(double * r, int * ind, double * coef, int * ibasis, int nterms, int npairs, int ncoefs, int isurf, double * v, double * dvdr)
{
  __shared__ double cache_v[_SIZE_BLOCK];
  __shared__ double cache_dvdr[_NUM_DF_PER * _SIZE_BLOCK];
 
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  int cache_id = threadIdx.x;

  int pid = blockIdx.y * blockDim.y + threadIdx.y;
  
  // useful work
  
  double v_temp = 0.0;
  double dvdr_temp[_NUM_DF_PER] = {0.0};

  while (id < nterms) {

    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow_int(exp(-r[j]), ind[id*npairs+j]);
    }
    arg *= coef[(isurf-1)*ncoefs + ibasis[id]-1];
    //printf("coef[%d] == %E\n", (isurf-1)*ncoefs + ibasis[id]-1, coef[(isurf-1)*ncoefs + ibasis[id]-1]);

    for(int j=0; j<_NUM_DF_PER; ++j) {
      const int indx_df = pid * _NUM_DF_PER + j;
      if(indx_df < npairs) dvdr_temp[j] -= arg*ind[id*npairs + indx_df];
    }
    
    v_temp += arg;
    
    id += blockDim.x * gridDim.x;
  }

  // set thread-local value
  
  cache_v[cache_id] = v_temp;
  for(int j=0; j<_NUM_DF_PER; ++j) {
    cache_dvdr[j*_SIZE_BLOCK + cache_id] = dvdr_temp[j];
  }

  // block
  
  __syncthreads();

  // manually reduce values from threads within block to master thread's value
  
  int i=blockDim.x / 2;
  while(i != 0) {
    if(cache_id < i) {
      cache_v[cache_id] += cache_v[cache_id + i];
      for(int j=0; j<_NUM_DF_PER; ++j) {
        cache_dvdr[j*_SIZE_BLOCK + cache_id] += cache_dvdr[j*_SIZE_BLOCK + cache_id + i];
      }
    }
    __syncthreads();
    i /= 2;
  }

  // store master thread's value in global array for host
  // energy & forces written blockDim.y times to same addresses
  
  if(cache_id == 0) {
    if(pid == 0) v[blockIdx.x] = cache_v[0];
    for(int j=0; j<_NUM_DF_PER; ++j) {
      dvdr[(pid*_NUM_DF_PER+j)*gridDim.x + blockIdx.x] = cache_dvdr[j*_SIZE_BLOCK];
    }
  }
}

__global__ void _compute_pip_algo2(double * r, int * ind, double * coef, int * ibasis, int nterms, int npairs, int ncoefs, int isurf, double * v)
{
  __shared__ double cache_v[_SIZE_BLOCK];
 
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  int cache_id = threadIdx.x;

  int pid = blockIdx.y * blockDim.y + threadIdx.y;
  
  // useful work
  
  double v_temp = 0.0;
  while (id < nterms) {

    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow_int(exp(-r[j]), ind[id*npairs+j]);
    }
    arg *= coef[(isurf-1)*ncoefs + ibasis[id]-1];

    v_temp += arg;
    
    id += blockDim.x * gridDim.x;
  }

  // set thread-local value
  
  cache_v[cache_id] = v_temp;

  // block
  
  __syncthreads();

  // manually reduce values from threads within block to master thread's value
  
  int i=blockDim.x / 2;
  while(i != 0) {
    if(cache_id < i) {
      cache_v[cache_id] += cache_v[cache_id + i];
    }
    __syncthreads();
    i /= 2;
  }

  // store master thread's value in global array for host
  // energy & forces written blockDim.y times to same addresses
  
  if(cache_id == 0) {
    if(pid == 0) v[blockIdx.x] = cache_v[0];
  }
}

#else

__global__ void _compute_pip(double * r, int * ind, double * coef, int * ibasis, int nterms, int npairs, int ncoefs, int isurf, double * v, double * dvdr)
{
  __shared__ double cache_v[_SIZE_BLOCK];
  __shared__ double cache_dvdr[_MAX_ATOM*(_MAX_ATOM-1)/2 * _SIZE_BLOCK];
 
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  int cache_id = threadIdx.x;

  // useful work
  
  double v_temp = 0.0;
  double dvdr_temp[_MAX_ATOM * (_MAX_ATOM-1) / 2] = {0.0};
  while (id < nterms) {

    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow_int(exp(-r[j]), ind[id*npairs+j]);
    }
    arg *= coef[(isurf-1)*ncoefs + ibasis[id]-1];
    
    for(int j=0; j<npairs; ++j) {
      dvdr_temp[j] -= arg*ind[id*npairs+j];
    }
    
    v_temp += arg;
    
    id += blockDim.x * gridDim.x;
  }

  // set thread-local value
  
  cache_v[cache_id] = v_temp;
  for(int j=0; j<npairs; ++j) {
    cache_dvdr[j*_SIZE_BLOCK + cache_id] = dvdr_temp[j];
  }

  // block
  
  __syncthreads();

  // manually reduce values from threads within block to master thread's value
  
  int i=blockDim.x / 2;
  while(i != 0) {
    if(cache_id < i) {
      cache_v[cache_id] += cache_v[cache_id + i];
      for(int j=0; j<npairs; ++j) {
        cache_dvdr[j*_SIZE_BLOCK + cache_id] += cache_dvdr[j*_SIZE_BLOCK + cache_id + i];
      }
    }
    __syncthreads();
    i /= 2;
  }

  // store master thread's value in global array for host
  
  if(cache_id == 0) {
    v[blockIdx.x] = cache_v[0];
    for(int j=0; j<npairs; ++j) {
      dvdr[j*blockDim.x + blockIdx.x] = cache_dvdr[j*_SIZE_BLOCK];
    }
  }
}

__global__ void _compute_pip(double * r, int * ind, double * coef, int * ibasis, int nterms, int npairs, int ncoefs, int isurf, double * v)
{
  __shared__ double cache_v[_SIZE_BLOCK];
 
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  int cache_id = threadIdx.x;

  // useful work
  
  double v_temp = 0.0;
  while (id < nterms) {

    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow_int(exp(-r[j]), ind[id*npairs+j]);
    }
    arg *= coef[(isurf-1)*ncoefs + ibasis[id]-1];
    
    v_temp += arg;
    
    id += blockDim.x * gridDim.x;
  }

  // set thread-local value
  
  cache_v[cache_id] = v_temp;

  // block
  
  __syncthreads();

  // manually reduce values from threads within block to master thread's value
  
  int i=blockDim.x / 2;
  while(i != 0) {
    if(cache_id < i) {
      cache_v[cache_id] += cache_v[cache_id + i];
    }
    __syncthreads();
    i /= 2;
  }

  // store master thread's value in global array for host
  
  if(cache_id == 0) {
    v[blockIdx.x] = cache_v[0];
  }
}

#endif

/* ---------------------------------------------------------------------- */

// maxpairs = maxatom * (maxatom-1)/2

// r(maxpairs)
// ind(mterm, maxpairs)
// coef(mcoef)
// ibasis(mterm)

// original indexing : 
   //      arg *= pow( exp(-r[j]), (double)ind[i+j*mterm] );

/* ---------------------------------------------------------------------- */

void Device::setup_device()
{
  grid_size = (nterms+_SIZE_BLOCK-1) / _SIZE_BLOCK;
  block_size = _SIZE_BLOCK;

#ifdef _ALGO2
  grid_y_size = (npairs+_NUM_DF_PER-1) / _NUM_DF_PER;
  block_y_size = 1;
  
  printf("Launching kernels w/ grid_size= (%lu, %lu) block_size= (%lu, %lu)\n",
	 grid_size, grid_y_size,
	 block_size, block_y_size);

  npmax = MAX(npairs, grid_y_size*block_y_size); // should just set to grid_y_size*block_y_size
#else
  grid_y_size = 1;
  block_y_size = 1;
  
  printf("Launching kernels w/ grid_size= (%lu, %lu) block_size= (%lu, %lu)\n",
	 grid_size, grid_y_size,
	 block_size, block_y_size);
  
  npmax = npairs;
#endif
  
  v_partial = (double*) dev_malloc_host(grid_size*sizeof(double));
  dvdr_partial = (double*) dev_malloc_host(npmax*grid_size*sizeof(double));
  
  d_v_partial = (double*) dev_malloc(grid_size*sizeof(double));
  d_dvdr_partial = (double*) dev_malloc(npmax*grid_size*sizeof(double));

  dev_stream_create(stream);
}

/* ---------------------------------------------------------------------- */

void Device::time_compute(double * r, double * vptr, double * dvdr, int isurf)
{ 
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  //  using std::chrono::milliseconds;
  using std::chrono::microseconds;
  
  // evaluate pip kernel
  
  for(int it=0; it<_NUM_ITER; ++it) {
    
    auto t1 = high_resolution_clock::now();

    compute(r, vptr, dvdr, isurf);
    
    // End the timer and print the time elapsed
    auto t2 = high_resolution_clock::now();
    //    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    auto ms_int = duration_cast<microseconds>(t2 - t1);
    int time = ms_int.count();
    
    printf(" C-Kernel : E= %f kcal/mol  time= %i (us)\n",*vptr, time);

    total_t += time;
  }
  
}

/* ---------------------------------------------------------------------- */

void Device::compute(double * r, double * vptr, double * dvdr, int isurf)
{ 
  // update coordinates on gpu

  dev_push(d_r, r, size_r);
  
  // evaluate pip kernel
    
#if defined(_CUDA_NVTX)
  nvtxRangePushA("compute");
#endif

#ifdef _ALGO2
  dim3 grid(grid_size, grid_y_size);
  dim3 block(block_size, block_y_size);
  
  _compute_pip_algo2<<<grid, block, 0, stream>>>(d_r,d_ind,d_coef,d_ibasis,nterms,npairs,ncoefs,isurf,d_v_partial,d_dvdr_partial);
#else
  _compute_pip<<<grid_size, block_size, 0, stream>>>(d_r,d_ind,d_coef,d_ibasis,nterms,npairs,ncoefs,isurf,d_v_partial,d_dvdr_partial);
#endif
  _CUDA_CHECK_ERRORS();
    
  dev_pull_async(d_v_partial, v_partial, grid_size*sizeof(double), stream);
  dev_pull_async(d_dvdr_partial, dvdr_partial, npmax*grid_size*sizeof(double), stream);

  dev_stream_wait(stream);

#if defined(_CUDA_NVTX)
  nvtxRangePop();  
  nvtxRangePushA("host_accumulate");
#endif

  double v = 0.0;
  for(int i=0; i<grid_size; ++i) v+= v_partial[i];

  *vptr = v;

  // use OpenMP here?
  // this might be large enough that we want to reduce on gpu
  
  for(int i=0; i<npairs; ++i) {
    double val = 0.0;
    for(int j=0; j<grid_size; ++j) val += dvdr_partial[i*grid_size + j];
    dvdr[i] = val;
  }
  
#if defined(_CUDA_NVTX)
  nvtxRangePop();
#endif
}

/* ---------------------------------------------------------------------- */

void Device::compute(double * r, double * vptr, int isurf)
{ 
  // update coordinates on gpu

  dev_push(d_r, r, size_r);
  
  // evaluate pip kernel
    
#if defined(_CUDA_NVTX)
  nvtxRangePushA("compute");
#endif

#ifdef _ALGO2
  dim3 grid(grid_size, grid_y_size);
  dim3 block(block_size, block_y_size);
  
  _compute_pip_algo2<<<grid, block, 0, stream>>>(d_r,d_ind,d_coef,d_ibasis,nterms,npairs,ncoefs,isurf,d_v_partial);
#else
  _compute_pip<<<grid_size, block_size, 0, stream>>>(d_r,d_ind,d_coef,d_ibasis,nterms,npairs,ncoefs,isurf,d_v_partial);
#endif
  _CUDA_CHECK_ERRORS();
    
  dev_pull_async(d_v_partial, v_partial, grid_size*sizeof(double), stream);

  dev_stream_wait(stream);

#if defined(_CUDA_NVTX)
  nvtxRangePop();  
  nvtxRangePushA("host_accumulate");
#endif

  double v = 0.0;
  for(int i=0; i<grid_size; ++i) v+= v_partial[i];

  *vptr = v;

#if defined(_CUDA_NVTX)
  nvtxRangePop();
#endif
}

#endif
