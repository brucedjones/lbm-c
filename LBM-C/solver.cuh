#ifndef SOLVER_H
#define SOLVER_H

#include "data_types.cuh"

// CUDA KERNEL PROTOTYPES
__global__ void iterate_bulk_kernel (Lattice *lattice, Domain *domain, bool store_macros);
__global__ void iterate_boundary_kernel (Lattice *lattice, Domain *domain, int offset, bool store_macros);
__global__ void iterate_all_kernel (Lattice *lattice, Domain *domain, int offset, int type);

// KERNEL SUPPORT DEVICE FUNCTION PROTOTYPES
__device__ inline int2 compute_boundary_coords(int idx, Domain *domain);

#endif