/* -*- c++ -*- */

#if defined(_USE_CPU)

#include "device.h"

#include <stdio.h>

/* ---------------------------------------------------------------------- */

inline double pow_int(double a, int n)
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
  num_threads = 1;
#pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }
  
  const int date = _OPENMP;

  double version;
  if     (date == 201107) version = 3.1;
  else if(date == 201307) version = 4.0;
  else if(date == 201511) version = 4.5;
  else if(date == 201611) version = 5.0;
  else if(date == 201811) version = 5.0;
  else if(date == 202011) version = 5.1;
  else {
    printf("[PM] Error: unknown omp version: omp_data= %i.\n",date);
    exit(1);
  }
  
  printf("\n  Using OPENMP v%3.1f\n",version);
  printf("  num_threads=     %i\n",num_threads);
  
  v_partial = (double*) dev_malloc_host(num_threads*sizeof(double));
  dvdr_partial = (double*) dev_malloc_host(npairs*num_threads*sizeof(double));
  
  //  d_v_partial = (double*) dev_malloc(grid_size*sizeof(double));
  //  d_dvdr_partial = (double*) dev_malloc(npairs*grid_size*sizeof(double));
}

/* ---------------------------------------------------------------------- */

void Device::time_compute(double * r, double * vptr, double * dvdr)
{ 
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  //  using std::chrono::milliseconds;
  using std::chrono::microseconds;

  // evaluate pip kernel
  
  for(int it=0; it<_NUM_ITER; ++it) {
    
    auto t1 = high_resolution_clock::now();

    compute(r, vptr, dvdr);
    
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

void Device::compute(double * r, double * vptr, double * dvdr)
{
  
  double v = 0.0;
  for(int i=0; i<npairs; ++i) dvdr[i] = 0.0;
  
  // evaluate pip kernel

#if 0  
  for(int i=0; i<nterms; ++i) {
    double arg = 1.0;
    for(int j=0; j<npairs; ++j) {
      arg *= pow_int( exp(-d_r[j]), d_ind[i*npairs+j] );
    }
    arg *= d_coef[d_ibasis[i]-1];

    for(int j=0; j<npairs; ++j) {
      dvdr[j] -= arg * d_ind[i*npairs+j];
    }
    
    v += arg;
  }

#else
  
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();

    int delta = nterms / num_threads;
    int rem = nterms - num_threads * delta;

    int from = (tid < rem) ? tid * (delta+1) : tid * delta + rem;
    int to = (tid < rem) ? from + delta + 1 : from + delta;

    double * dvdr_ = &(dvdr_partial[tid * npairs]);

    double v_ = 0.0;
    for(int j=0; j<npairs; ++j) dvdr_[j] = 0.0;

    for(int i=from; i<to; ++i) {
      double arg = 1.0;
      for(int j=0; j<npairs; ++j) arg *= pow_int( exp(-d_r[j]), d_ind[i*npairs+j] );
      arg *= d_coef[d_ibasis[i]-1];
      
      for(int j=0; j<npairs; ++j) dvdr_[j] -= arg * d_ind[i*npairs+j];
      
      v_ += arg;
    }

#pragma omp atomic
    v += v_;
    
#pragma omp barrier
    
    // Reduce thread-local values

    delta = npairs / num_threads;
    rem = npairs - num_threads * delta;

    from = (tid < rem) ? tid * (delta+1) : tid * delta + rem;
    to = (tid < rem) ? from + delta + 1 : from + delta;

    for(int i=from; i<to; ++i) {
      for(int j=0; j<num_threads; ++j) dvdr[i] += dvdr_partial[j*npairs+i];
    }
  }
  
#endif
  
  *vptr = v;
}

#endif
