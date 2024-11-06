#if defined(_GPU_OPENMP)

#include <stdio.h>
#include <iostream>

#include "pm_openmp.h"

int dev_num_devices()
{
  int num_devices = omp_get_num_devices();
  _OMP_CHECK_ERRORS();
  
  return num_devices;
}

void dev_properties(int ndev)
{
  int num_devices = omp_get_num_devices();
  int default_device = omp_get_default_device();
  int host = omp_get_initial_device();

  int num_teams = -1;
  int num_threads = -1;
#pragma omp target teams map(tofrom: num_teams, num_threads)
  {
    num_teams = omp_get_num_teams();
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
  printf("  num_devices=     %i\n",num_devices);
  printf("  Default device=  %i\n",default_device);
  printf("  Host=            %i\n",host);
  printf("  num_teams=       %i\n",num_teams);
  printf("  num_threads=     %i\n",num_threads);
}

int dev_check_peer(int rank, int ngpus)
{
  int err = 0;
  // if(rank == 0) printf("\nChecking P2P Access\n");
  // for(int ig=0; ig<ngpus; ++ig) {
  //   cudaSetDevice(ig);
  //   //if(rank == 0) printf("Device i= %i\n",ig);

  //   int n = 1;
  //   for(int jg=0; jg<ngpus; ++jg) {
  //     if(jg != ig) {
  //       int access;
  //       cudaDeviceCanAccessPeer(&access, ig, jg);
  //       n += access;

  //       //if(rank == 0) printf("  --  Device j= %i  access= %i\n",jg,access);
  //     }
  //   }
  //   if(n != ngpus) err += 1;
  // }

  return err;
}

void dev_set_device(int id)
{
  omp_set_default_device(id);
  _OMP_CHECK_ERRORS();
}

void * dev_malloc(int N)
{
  int id = omp_get_default_device();
  void * ptr = omp_target_alloc(N, id);
  _OMP_CHECK_ERRORS();
  return ptr;
}

void * dev_malloc_host(int N)
{
  void * ptr = omp_alloc(N, omp_default_mem_alloc);
  _OMP_CHECK_ERRORS();
  return ptr;
}

void dev_free(void * ptr)
{
  int id = omp_get_default_device();
  omp_target_free(ptr, id);
  _OMP_CHECK_ERRORS();
}

void dev_free_host(void * ptr)
{
  omp_free(ptr, omp_default_mem_alloc);
  _OMP_CHECK_ERRORS();
}

void dev_push(void * d_ptr, void * h_ptr, int N)
{
  int gpu = omp_get_default_device();
  int host = omp_get_initial_device();
  omp_target_memcpy(d_ptr, h_ptr, N, 0, 0, gpu, host);
  _OMP_CHECK_ERRORS();
}

void dev_pull(void * d_ptr, void * h_ptr, int N)
{
  int gpu = omp_get_default_device();
  int host = omp_get_initial_device();
  omp_target_memcpy(h_ptr, d_ptr, N, 0, 0, host, gpu);
  _OMP_CHECK_ERRORS();
}

void dev_copy(void * a_ptr, void * b_ptr, int N)
{
  int gpu = omp_get_default_device();
  omp_target_memcpy(b_ptr, a_ptr, N, 0, 0, gpu, gpu);
  _OMP_CHECK_ERRORS();
}
#endif
