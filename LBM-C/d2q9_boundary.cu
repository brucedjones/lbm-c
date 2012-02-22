#ifndef D2Q9_BOUNDARY
#define D2Q9_BOUNDARY

// Necessary includes
#include "macros.cu"
#include "d2q9_boundary.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ __noinline__ void zh_pressure_x(Node *current_node, double *rho_boundary)
{
	// COMPUTE MACROS
	current_node->rho = *rho_boundary;
	current_node->u[0] = 1.0-((1.0/current_node->rho)*(current_node->f[0]+current_node->f[2]+current_node->f[4]+2.0*(current_node->f[3]+current_node->f[6]+current_node->f[7])));
	current_node->u[1] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[1] = current_node->f[3] + ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[5] = current_node->f[7] - ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) + ((1.0/6.0)*current_node->rho*current_node->u[0]);
	current_node->f[8] = current_node->f[6] + ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) + ((1.0/6.0)*current_node->rho*current_node->u[0]);

}

__device__ __noinline__ void zh_pressure_X(Node *current_node, double *rho_boundary)
{
	// COMPUTE MACROS
	current_node->rho = *rho_boundary;
	current_node->u[0] = -1.0+((1.0/current_node->rho)*(current_node->f[0]+current_node->f[2]+current_node->f[4]+2.0*(current_node->f[1]+current_node->f[5]+current_node->f[8])));
	current_node->u[1] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[3] = current_node->f[1] - ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[6] = current_node->f[8] - ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) - ((1.0/6.0)*current_node->rho*current_node->u[0]);
	current_node->f[7] = current_node->f[5] + ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) - ((1.0/6.0)*current_node->rho*current_node->u[0]);
}

#endif
