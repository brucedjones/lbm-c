#ifndef COLLISION_H
#define COLLISION_H

// SUPPORT FUNCTIONS
__device__ inline double u_square(Node *current_node);
__device__ inline double e_mul_u(Node *current_node, int **e, int *i);
__device__ inline void turbulent_viscosity(Node *current_node, double *f_eq, double *tau);

// BOUNDARY CONDITION DEVICE FUNCTION PROTOTYPES
__device__ void bgk_collision(Node *current_node, double *tau);
__device__ void bgk_guo_collision(Node *current_node, double *tau);
__device__ void bgk_ntpor_collision(Node *current_node, double *tau);
__device__ void bgk_ntpor_guo_collision(Node *current_node, double *tau);
__device__ void mrt_collision(Node *current_node, double *tau);
__device__ void mrt_guo_collision(Node *current_node, double *tau);
__device__ void mrt_ntpor_collision(Node *current_node, double *tau);
__device__ void mrt_ntpor_guo_collision(Node *current_node, double *tau);
__device__ void meq_d2q9(Node *current_node, double *meq);
__device__ void meq_d3q15(Node *current_node, double *meq);
__device__ void bounceback(Node *current_node, double *tau);

#endif