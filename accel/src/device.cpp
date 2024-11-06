/* -*- c++ -*- */

#include <stdio.h>
#include "device.h"

#if defined(_CUDA_NVTX)
#include "nvToolsExt.h"
#endif

/* ---------------------------------------------------------------------- */

Device::Device()
{
  printf("LIBGPU: created device\n");
  total_t = 0.0;

  d_r = nullptr;
  d_ind = nullptr;
  d_ibasis = nullptr;
  d_coef = nullptr;

  v_partial = nullptr;
  dvdr_partial = nullptr;
  d_v_partial = nullptr;
  d_dvdr_partial = nullptr;
}

/* ---------------------------------------------------------------------- */

Device::~Device()
{
  printf("LIBGPU: destroying device\n");
  
#if defined(_CUDA_NVTX)
  nvtxRangePushA("cleanup");
#endif
  
  dev_free(d_r);
  dev_free(d_ind);
  dev_free(d_ibasis);
  dev_free(d_coef);
  
  dev_free(d_v_partial);
  dev_free(d_dvdr_partial);

  dev_free_host(v_partial);
  dev_free_host(dvdr_partial);
  
#if defined(_CUDA_NVTX)
  nvtxRangePop();
#endif
}

/* ---------------------------------------------------------------------- */

int Device::get_num_devices()
{
  printf("LIBGPU: getting number of devices\n");
  return dev_num_devices();
}

/* ---------------------------------------------------------------------- */
    
void Device::set_device(int id)
{
  printf("LIBGPU: setting device id= %i\n",id);
  dev_set_device(id);
}

/* ---------------------------------------------------------------------- */

double Device::host_compute_pip(double * r, int * ind, double * coef, int * ibasis, int nterms, int npairs)
{
  
  double v = 0.0;
  for(int i=0; i<nterms; ++i) {
    
    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow( exp(-r[j]), ind[i*npairs+j] );
    }
    v += arg * coef[ibasis[i]-1];
  }

  return v;
}

/* ---------------------------------------------------------------------- */

void Device::hist_exponent(int * ind, int nterms, int npairs)
{
  int max = 10;

  int * hist = (int*) malloc(max*sizeof(int));

  for(int i=0; i<max; ++i) hist[i] = 0;
  
  for(int i=0; i<nterms; ++i)
    for(int j=0; j<npairs; ++j) {
      int k = ind[i*npairs+j];
      if(k >= max) {
	printf("ind >= 10 detected: ij= %i %i  ind= %i\n",i,j,k);
	exit(1);
      } else if(k < 0) {
	printf("ind < 0 detected: ij= %i %i  ind= %i\n",i,j,k);
      } else hist[ k ]++;
    }

  printf("\nHistogram of exponents: nterms= %i  npairs= %i\n",nterms,npairs);
  for(int i=0; i<max; ++i) printf("i= %i  hist= %i\n",i,hist[i]);

  double frac = 0.0;
  for(int i=0; i<5; ++i) frac+= hist[i];
  frac /= nterms*npairs;
  printf("%% of hits to first 5 values = %f %%\n\n",frac*100);
  
  free(hist);
}

/* ---------------------------------------------------------------------- */

void Device::setup(int * ind, double * coef, int * ibasis, int _nterms, int _ncoefs, int _npairs, int _nsurfs)
{
  // save local copies

  nterms = _nterms;
  ncoefs = _ncoefs;
  npairs = _npairs;
  nsurfs = _nsurfs;
  
  // create device buffers

  size_r = npairs*sizeof(double);
  size_ind = npairs*nterms*sizeof(int);
  size_ibasis = nterms*sizeof(int);
  size_coef = nsurfs*ncoefs*sizeof(double);

#if defined(_CUDA_NVTX)
  nvtxRangePushA("setup");
#endif
  
  d_r = (double*) dev_malloc(size_r);
  d_ind = (int*) dev_malloc(size_ind);
  d_ibasis = (int*) dev_malloc(size_ibasis);
  d_coef = (double*) dev_malloc(size_coef);

  dev_push(d_ind, ind, size_ind);
  dev_push(d_ibasis, ibasis, size_ibasis);
  dev_push(d_coef, coef, size_coef);
  
  // histogram of integer exponents

  // hist_exponent(ind, nterms, npairs);

  // device-specific setup

  setup_device();
  
#if defined(_CUDA_NVTX)
  nvtxRangePop();
#endif
  
  //double _v = host_compute_pip(r, ind, coef, ibasis, nterms, npairs);
  //printf(" C-Kernel Reference E= %f\n\n",_v); 
}

/* ---------------------------------------------------------------------- */

