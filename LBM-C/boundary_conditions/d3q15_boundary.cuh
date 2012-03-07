#ifndef D3Q15_BOUNDARY_H
#define D3Q15_BOUNDARY_H

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ __noinline__ void zh_pressure_x(Node *current_node, double *rho_boundary);
__device__ __noinline__ void zh_pressure_X(Node *current_node, double *rho_boundary);

#endif