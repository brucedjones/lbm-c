#ifndef COLLISION
#define COLLISION

// Necessary includes
#include "macros.cu"
#include "collision.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ __noinline__ void bgk_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B)
{
	double f_eq, u_sq, cu;

	u_sq = 1.5*(current_node->ux*current_node->ux + current_node->uy*current_node->uy);
	for(int i=0;i<Q;i++)
	{
		cu = 3.0*(ex[i]*current_node->ux+ey[i]*current_node->uy);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/tau) * (current_node->f[i]-f_eq);
	}
}

__device__ __noinline__ void guo_bgk_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B)
{
	double f_eq, u_sq, cu, F_coeff[DIM], force_term;
	
	current_node->ux = current_node->ux + (1/(2*current_node->rho))*current_node->F[0];
	current_node->uy = current_node->uy + (1/(2*current_node->rho))*current_node->F[1];

	u_sq = 1.5*(current_node->ux*current_node->ux + current_node->uy*current_node->uy);

	for(int i=0;i<Q;i++)
	{
		F_coeff[0] = omega[i]*(1-(1/(2*tau)))*(((ex[i]-current_node->ux)*3)+(ex[i]*9*((ex[i]*current_node->ux)+(ey[i]*current_node->uy))));
		F_coeff[1] = omega[i]*(1-(1/(2*tau)))*(((ey[i]-current_node->uy)*3)+(ey[i]*9*((ex[i]*current_node->ux)+(ey[i]*current_node->uy))));

		force_term = F_coeff[0]*current_node->F[0]+F_coeff[1]*current_node->F[1];

		cu = 3.0*(ex[i]*current_node->ux+ey[i]*current_node->uy);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/tau) * (current_node->f[i]-f_eq)+force_term;
	}
}

__device__ __noinline__ void nt_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B)
{
	double f_eq, u_sq, cu, collision_bgk, collision_s;

	u_sq = 1.5*(current_node->ux*current_node->ux + current_node->uy*current_node->uy);
	for(int i=0;i<Q;i++)
	{
		cu = 3.0*(ex[i]*current_node->ux+ey[i]*current_node->uy);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		collision_bgk = (1.0/tau) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];

		current_node->f[i] = current_node->f[i] - (1-B)*collision_bgk + B*collision_s;
	}
}

__device__ void guo_nt_collision(Node *current_node, int opp[Q], int ex[Q], int ey[Q], double omega[Q], double tau, double B)
{
	double f_eq, u_sq, cu, collision_bgk, collision_s, F_coeff[DIM], force_term;

	current_node->ux = current_node->ux + (1/(2*current_node->rho))*current_node->F[0];
	current_node->uy = current_node->uy + (1/(2*current_node->rho))*current_node->F[1];

	u_sq = 1.5*(current_node->ux*current_node->ux + current_node->uy*current_node->uy);

	for(int i=0;i<Q;i++)
	{
		F_coeff[0] = omega[i]*(1-(1/(2*tau)))*(((ex[i]-current_node->ux)*3)+(ex[i]*9*((ex[i]*current_node->ux)+(ey[i]*current_node->uy))));
		F_coeff[1] = omega[i]*(1-(1/(2*tau)))*(((ey[i]-current_node->uy)*3)+(ey[i]*9*((ex[i]*current_node->ux)+(ey[i]*current_node->uy))));

		force_term = F_coeff[0]*current_node->F[0]+F_coeff[1]*current_node->F[1];

		cu = 3.0*(ex[i]*current_node->ux+ey[i]*current_node->uy);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		collision_bgk = (1.0/tau) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];

		current_node->f[i] = current_node->f[i] - (1-B)*(collision_bgk+force_term) + B*collision_s;
	}
}

#endif
