#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#ifdef _OPENMP
#include <pthread.h>
#include <omp.h>

#else

static inline int omp_get_num_procs(void){
  return 1;
}

static inline int omp_get_thread_limit(void){
  return 1;
}

static inline int omp_get_max_threads(void){
  return 1;
}

static inline int omp_get_thread_num(void) {
  return 0;
}

static inline int omp_in_parallel(void) {
  return 0;
}

#endif
