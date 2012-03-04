#ifndef CUDA_UTIL_H
#define CUDA_UTIL_H

template<class T>
void combi_malloc(T**, T**, size_t);
void cudasafe( cudaError_t, char*);
void Check_CUDA_Error(const char*);

#endif