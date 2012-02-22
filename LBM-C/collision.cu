#ifndef COLLISION
#define COLLISION

// Necessary includes
#include "macros.cu"
#include "collision.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ __noinline__ void bgk_collision(Node *current_node, int *opp, int *ex, int *ey, double *omega, double *tau, double *B)
{
	double f_eq, u_sq, cu;

	u_sq = 1.5*(current_node->u[0]*current_node->u[0] + current_node->u[1]*current_node->u[1]);
	for(int i=0;i<Q;i++)
	{
		cu = 3.0*(ex[i]*current_node->u[0]+ey[i]*current_node->u[1]);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq);
	}
}

__device__ __noinline__ void guo_bgk_collision(Node *current_node, int *opp, int *ex, int *ey, double *omega, double *tau, double *B)
{
	double f_eq, u_sq, cu, F_coeff[DIM], force_term;
	
	current_node->u[0] = current_node->u[0] + (1/(2*current_node->rho))*current_node->F[0];
	current_node->u[1] = current_node->u[1] + (1/(2*current_node->rho))*current_node->F[1];

	u_sq = 1.5*(current_node->u[0]*current_node->u[0] + current_node->u[1]*current_node->u[1]);

	for(int i=0;i<Q;i++)
	{
		F_coeff[0] = omega[i]*(1-(1/(2*(*tau))))*(((ex[i]-current_node->u[0])*3)+(ex[i]*9*((ex[i]*current_node->u[0])+(ey[i]*current_node->u[1]))));
		F_coeff[1] = omega[i]*(1-(1/(2*(*tau))))*(((ey[i]-current_node->u[1])*3)+(ey[i]*9*((ex[i]*current_node->u[0])+(ey[i]*current_node->u[1]))));

		force_term = F_coeff[0]*current_node->F[0]+F_coeff[1]*current_node->F[1];

		cu = 3.0*(ex[i]*current_node->u[0]+ey[i]*current_node->u[1]);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq)+force_term;
	}
}

__device__ __noinline__ void ntpor_collision(Node *current_node, int *opp, int *ex, int *ey, double *omega, double *tau, double *B)
{
	double f_eq, u_sq, cu, collision_bgk, collision_s, tmp[Q];

	u_sq = 1.5*(current_node->u[0]*current_node->u[0] + current_node->u[1]*current_node->u[1]);
	for(int i=0;i<Q;i++)
	{
		cu = 3.0*(ex[i]*current_node->u[0]+ey[i]*current_node->u[1]);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];
		
		tmp[i] = current_node->f[i] - (1-(*B))*collision_bgk + (*B)*collision_s;
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}

}

__device__ void guo_ntpor_collision(Node *current_node, int *opp, int *ex, int *ey, double *omega, double *tau, double *B)
{
	double f_eq, u_sq, cu, collision_bgk, collision_s, F_coeff[DIM], force_term, tmp[Q];

	current_node->u[0] = current_node->u[0] + (1/(2*current_node->rho))*current_node->F[0];
	current_node->u[1] = current_node->u[1] + (1/(2*current_node->rho))*current_node->F[1];

	u_sq = 1.5*(current_node->u[0]*current_node->u[0] + current_node->u[1]*current_node->u[1]);

	for(int i=0;i<Q;i++)
	{
		F_coeff[0] = omega[i]*(1-(1/(2*(*tau))))*(((ex[i]-current_node->u[0])*3)+(ex[i]*9*((ex[i]*current_node->u[0])+(ey[i]*current_node->u[1]))));
		F_coeff[1] = omega[i]*(1-(1/(2*(*tau))))*(((ey[i]-current_node->u[1])*3)+(ey[i]*9*((ex[i]*current_node->u[0])+(ey[i]*current_node->u[1]))));

		force_term = F_coeff[0]*current_node->F[0]+F_coeff[1]*current_node->F[1];

		cu = 3.0*(ex[i]*current_node->u[0]+ey[i]*current_node->u[1]);
		f_eq = current_node->rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];

		tmp[i] = current_node->f[i] - (1-(*B))*(collision_bgk+force_term) + (*B)*collision_s;
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}
}

__device__ void bounceback(Node *current_node, int *opp, int *ex, int *ey, double *omega, double *tau, double *B)
{
	double tmp[Q];
	for(int i=0;i<Q;i++)
	{
		tmp[i] = current_node->f[i];
	}

	for(int i=0;i<Q;i++)
	{
		current_node->f[i] = tmp[opp[i]];
	}

	current_node->u[0] = 0;
	current_node->u[1] = 0;
	current_node->rho = 0;
}

#endif
