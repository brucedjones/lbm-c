#ifndef SOLVER_H
#define SOLVER_H

#include "data_types.cuh"

// CUDA KERNEL PROTOTYPES
__global__ void iterate_kernel (Lattice *lattice, Domain *domain_arrays, bool store_macros, int t);

#endif