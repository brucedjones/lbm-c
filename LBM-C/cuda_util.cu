#ifndef CUDA_UTIL
#define CUDA_UTIL

#include "cuda_util.cuh"

template<class T>
void combi_malloc(T **host_pointer, T **device_pointer, size_t size)
{
	*host_pointer = (T *)malloc(size);
	cudasafe(cudaMalloc((void **)&*device_pointer,size), "Model Builder: Device memory allocation failed!");
}

// ERROR CHECKING FOR MEMORY ALLOCATION
void cudasafe( cudaError_t error, char* message)
{
   if(error!=cudaSuccess) { fprintf(stderr,"ERROR: %s : %s\n",message,cudaGetErrorString(error)); exit(-1); }
}

// ERROR CHECKING FOR KERNEL EXECUTION
void Check_CUDA_Error(const char *message)
{
   cudaError_t error = cudaGetLastError();
   if(error!=cudaSuccess) {
      fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
      exit(-1);
   }                         
}
#endif