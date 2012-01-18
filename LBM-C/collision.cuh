#ifndef D2Q9_BOUNDARY_H
#define D2Q9_BOUNDARY_H

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ Node bgk_collision(Node input, double rho_boundary);
__device__ Node nt_collision(Node input, double rho_boundary);

#endif