#ifndef SOLVER_H
#define SOLVER_H

#include "data_types.cuh"

// CUDA KERNEL PROTOTYPES
__global__ void iterate_bulk_kernel (Lattice *lattice_1, Lattice *lattice_2, Domain *domain);
__global__ void iterate_boundary_kernel (Lattice *lattice_1, Lattice *lattice_2, Domain *domain, int offset);
__global__ void iterate_all_kernel (Lattice *lattice_1, Lattice *lattice_2, Domain *domain, int offset, int type);

// KERNEL SUPPORT DEVICE FUNCTION PROTOTYPES
__device__ inline int2 compute_boundary_coords(int idx, Domain *domain);

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ Node zh_pressure_x(Node input, float rho_boundary);
__device__ Node zh_pressure_X(Node input, float rho_boundary);

#endif