#ifndef D2Q9_SF_DEFS_H
#define D2Q9_SF_DEFS_H

#include "d2q9_sf_defs.cu"

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ __noinline__ void sf_x(Node *current_node, Lattice *lattice);
__device__ __noinline__ void sf_X(Node *current_node, Lattice *lattice);
__device__ __noinline__ void sf_y(Node *current_node, Lattice *lattice);
__device__ __noinline__ void sf_Y(Node *current_node, Lattice *lattice);

#endif