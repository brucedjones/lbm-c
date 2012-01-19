#ifndef COLLISION_H
#define COLLISION_H

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ void bgk_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B);
__device__ void guo_bgk_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B);
__device__ void nt_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B);
__device__ void guo_nt_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B);

#endif