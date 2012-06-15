#ifndef MACRO_BC_H
#define MACRO_BC_H

__device__ __noinline__ void presc_ux(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_uy(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_uz(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_uy_ux(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_uz_ux(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_uz_uy(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_ux_uy_uz(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_ux(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_uy(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_uz(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_ux_uy(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_ux_uz(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_uy_uz(Node *current_node, Domain *domain);
__device__ __noinline__ void presc_rho_ux_uy_uz(Node *current_node, Domain *domain);

#endif