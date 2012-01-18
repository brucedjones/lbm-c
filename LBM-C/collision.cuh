#ifndef COLLISION_H
#define COLLISION_H

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ void bgk_collision(Node *current_node, int *opp, int *omega, int *ex, int *ey, double *F, double tau, double B);
__device__ void nt_collision(Node *current_node, int *opp, int *omega, int *ex, int *ey, double *F, double tau, double B);
__device__ void guo_nt_collision(Node *current_node, int *opp, int *omega, int *ex, int *ey, double *F, double tau, double B);

#endif