#ifndef D2Q9_ZH_DEFS_H
#define D2Q9_ZH_DEFS_H

#include "d2q9_zh_defs.cu"

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ __noinline__ void zh_pressure_x(Node *current_node, Lattice *lattice);
__device__ __noinline__ void zh_pressure_X(Node *current_node, Lattice *lattice);

#endif