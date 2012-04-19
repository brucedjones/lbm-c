#ifndef COLLISION
#define COLLISION

// Necessary includes
#include "macros.cu"
#include "collision.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ collision collision_functions[5] = { bgk_collision, guo_bgk_collision, ntpor_collision, guo_ntpor_collision, bounceback};

__device__ inline double u_square(Node *current_node)
{
	double value = 0;

	#pragma unroll
	for(int d = 0; d<DIM; d++)
	{
		value += (current_node->u[d]*current_node->u[d]);
	}

	return value*1.5;
}

__device__ inline double e_mul_u(Node *current_node, int e[DIM][Q], int *i)
{
	double value = 0;

	#pragma unroll
	for(int d = 0; d<DIM; d++)
	{
		value += (e[d][*i]*current_node->u[d]);
	}

	return value*3.;
}

__device__ __noinline__ void bgk_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B)
{
	double f_eq, u_sq, eu;

	u_sq = u_square(current_node);
	for(int i=0;i<Q;i++)
	{
		eu = e_mul_u(current_node, e, &i);
		f_eq = current_node->rho*omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq);
	}
}

__device__ __noinline__ void guo_bgk_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B)
{
	double f_eq, u_sq, eu, F_coeff[DIM], force_term;
	int d;
	
	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/2)*current_node->rho*current_node->F[d];
	}

	u_sq = u_square(current_node);

	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = omega[i]*(1-(1/(2*(*tau))))*(((e[d][i]-current_node->u[d])*3)+(e[d][i]*9*((e[0][i]*current_node->u[0])+(e[1][i]*current_node->u[1])+(e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = omega[i]*(1-(1/(2*(*tau))))*(((e[d][i]-current_node->u[d])*3)+(e[d][i]*9*((e[0][i]*current_node->u[0])+(e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term += F_coeff[d]*current_node->F[d];
		}

		eu = e_mul_u(current_node, e, &i);
		f_eq = current_node->rho*omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	
		current_node->f[i] = current_node->f[i] - (1.0/(*tau)) * (current_node->f[i]-f_eq)+force_term;
	}
}

__device__ __noinline__ void ntpor_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B)
{
	double f_eq, u_sq, eu, collision_bgk, collision_s, tmp[Q];

	u_sq = u_square(current_node);
	for(int i=0;i<Q;i++)
	{
		eu = e_mul_u(current_node, e, &i);
		f_eq = current_node->rho*omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];
		
		tmp[i] = current_node->f[i] - (1-(*B))*collision_bgk + (*B)*collision_s;
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}

}

__device__ void guo_ntpor_collision(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B)
{
	double f_eq, u_sq, eu, collision_bgk, collision_s, F_coeff[DIM], force_term, tmp[Q];
	int d;

	#pragma unroll
	for(d = 0; d<DIM; d++)
	{
		current_node->u[d] = current_node->u[d] + (1/2)*current_node->rho*current_node->F[d];
	}

	u_sq = u_square(current_node);

	for(int i=0;i<Q;i++)
	{
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
		#if DIM > 2
			F_coeff[d] = omega[i]*(1-(1/(2*(*tau))))*(((e[d][i]-current_node->u[d])*3)+(e[d][i]*9*((e[0][i]*current_node->u[0])+(e[1][i]*current_node->u[1])+(e[2][i]*current_node->u[2]))));
		#else
			F_coeff[d] = omega[i]*(1-(1/(2*(*tau))))*(((e[d][i]-current_node->u[d])*3)+(e[d][i]*9*((e[0][i]*current_node->u[0])+(e[1][i]*current_node->u[1]))));
		#endif
		}

		force_term = 0;
		#pragma unroll
		for(d = 0; d<DIM; d++)
		{
			force_term += F_coeff[d]*current_node->F[d];
		}

		eu = e_mul_u(current_node, e, &i);
		f_eq = current_node->rho*omega[i]*(1.0+eu+(0.5*eu*eu)-u_sq);
	
		collision_bgk = (1.0/(*tau)) * (current_node->f[i]-f_eq);
		collision_s = current_node->f[opp[i]]-current_node->f[i];

		tmp[i] = current_node->f[i] - (1-(*B))*(collision_bgk) + (*B)*collision_s + (1-(*B))*force_term;
	}

	for(int i =0;i<Q;i++)
	{
		current_node->f[i] = tmp[i];
	}
}

__device__ void bounceback(Node *current_node, int *opp, int e[DIM][Q], double *omega, double *tau, double *B)
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
