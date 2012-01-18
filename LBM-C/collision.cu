#ifndef D2Q9_BOUNDARY
#define D2Q9_BOUNDARY

// Necessary includes
#include "macros.cu"
#include "d2q9_boundary.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ __noinline__ Node bgk_collision(Node current_node, double rho_boundary)
{
	Node output; //= input;

	u_sq = 1.5*(current_node.ux*current_node.ux + current_node.uy*current_node.uy);
	for(int i=0;i<Q;i++)
	{
		cu = 3.0*(ex[i]*current_node.ux+ey[i]*current_node.uy);
		f_eq = current_node.rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
		current_node.f[i] = current_node.f[i] - (1.0/tau) * (current_node.f[i]-f_eq);
	}
	
	return output;
}

__device__ __noinline__ Node nt_collision(Node input, double rho_boundary)
{
	Node output; //= input;

	// COPY KNOWN f's
	output.f[0] = input.f[0];
	output.f[2] = input.f[2];
	output.f[4] = input.f[4];
	output.f[1] = input.f[1];
	output.f[5] = input.f[5];
	output.f[8] = input.f[8];

	// COMPUTE MACROS
	output.rho = rho_boundary;
	output.ux = -1.0+((1.0/output.rho)*(input.f[0]+input.f[2]+input.f[4]+2.0*(input.f[1]+input.f[5]+input.f[8])));
	output.uy = 0.0;

	// COMPUTE UNKNOWN f's
	output.f[3] = input.f[1] - ((2.0/3.0)*output.rho*output.ux);
	output.f[6] = input.f[8] - ((1.0/2.0)*(input.f[2]-input.f[4])) - ((1.0/6.0)*output.rho*output.ux);
	output.f[7] = input.f[5] + ((1.0/2.0)*(input.f[2]-input.f[4])) - ((1.0/6.0)*output.rho*output.ux);
	
	return output;
}

#endif
