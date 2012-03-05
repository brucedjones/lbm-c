#ifndef COLLISION_H
#define COLLISION_H

// SUPPORT FUNCTIONS
__device__ inline double u_square(Node *current_node);
__device__ inline double e_mul_u(Node *current_node, int **e, int *i);

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ void bgk_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B);
__device__ void guo_bgk_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B);
__device__ void ntpor_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B);
__device__ void guo_ntpor_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B);
__device__ void bounceback(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B);

#endif