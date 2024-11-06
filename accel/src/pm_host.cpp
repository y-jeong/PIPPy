#if defined(_USE_CPU)

#include <stdio.h>
#include <iostream>
#include <cstring>

#include <iomanip>
#include <vector>
#include <tuple>

#include "pm_host.h"

int dev_num_devices() {return 0;}

void dev_properties(int ndev) {}

int dev_check_peer(int rank, int ngpus) {return 0;}

void dev_set_device(int id) {}

int dev_get_device() {return 0;}

void * dev_malloc(size_t N) {return malloc(N);}

void * dev_malloc_host(size_t N) {return malloc(N);}

void dev_free(void * ptr) {free(ptr);}

void dev_free_host(void * ptr) {free(ptr);}

void dev_push(void * d_ptr, void * h_ptr, size_t N) {memcpy(d_ptr, h_ptr, N);}

void dev_pull(void * d_ptr, void * h_ptr, size_t N) {memcpy(h_ptr, d_ptr, N);}

void dev_copy(void * dest, void * src, size_t N) {memcpy(dest, src, N);}

void dev_check_pointer(int rnk, const char * name, void * ptr)
{
  if(ptr != nullptr) printf("(%i) ptr %s is hostPointer\n",rnk,name);
}

#endif
