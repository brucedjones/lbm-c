#ifndef SOLVER
#define SOLVER

#include "solver.cuh"
#include "d2q9_boundary.cu"
#include "collision.cu"

// LIST OF AVAILABLE BOUNDRY TYPES AND COLLISION FUNCTIONS
__device__ boundary_condition boundary_conditions[2] = { zh_pressure_x, zh_pressure_X};
__device__ collision collision_functions[5] = { bgk_collision, guo_bgk_collision, ntpor_collision, guo_ntpor_collision, bounceback};

__global__ void iterate_kernel (Lattice *lattice, DomainArray *domain_arrays, DomainConstant *domain_constants, bool store_macros)
{
	// Declare Variables
	double omega[Q], B;
	int i2d, ex[Q], ey[Q], opp[Q];
	int2 length;
	Node current_node;

	// Initialise variables
	LOAD_EX(ex);
	LOAD_EY(ey);
	LOAD_OMEGA(omega);
	LOAD_OPP(opp);
	current_node.rho = 0; current_node.ux = 0; current_node.uy = 0;
	
	// Compute coordinates
	int x = (blockDim.x*blockIdx.x)+threadIdx.x;
	int y = (blockDim.y*blockIdx.y)+threadIdx.y;

	// Load domain configuration
	length.x = domain_constants->length.x;
	length.y = domain_constants->length.y;
	int domain_size = length.x*length.y;
	double tau = domain_constants->tau;
	int i2d_prime = x + y*length.x;
	
	if(x<length.x && y<length.y)
	{
		// Set collision type and optional forces
		// The type specified in domain_constants must be multiplied by two to match the listing
		// order in the collision_functions array, an additional 1 is added to the collision type
		// to specify a collision with guo body forces
		int collision_modifier = 0;
		if(domain_constants->forcing==true)
		{
			current_node.F[0] = domain_arrays->force[i2d_prime];
			current_node.F[1] = domain_arrays->force[domain_size+i2d_prime];
			if(current_node.F[0] > 0 && current_node.F[1] > 0) collision_modifier = 1;
		}
		int collision_type = (domain_constants->collision_type*2)+collision_modifier;

		// Load boundary condition
		int boundary_type = domain_arrays->boundary_type[i2d_prime];
		double boundary_value = domain_arrays->boundary_value[i2d_prime];
	
		// Load Geometry
		B = domain_arrays->geometry[i2d_prime];
		if(B==1) collision_type = 4;
	
		// STREAMING - UNCOALESCED READ
		int target_x, target_y;
		for(int i = 0; i<Q; i++)
		{
			target_x = x+ex[i]; target_y = y+ey[i];
			//PERIODIC BOUNDARY
			if(target_x>(length.x-1)) target_x = 0; if(target_x<0) target_x = length.x-1;
			if(target_y>(length.y-1)) target_y = 0; if(target_y<0) target_y = length.y-1;
	
			i2d = (target_x + target_y*length.x)+opp[i]*(domain_size);
			
			// UNCOALESCED READ
			current_node.f[opp[i]] = lattice->f_prev[i2d];
	
			current_node.rho += current_node.f[opp[i]];
			current_node.ux += ex[opp[i]]*current_node.f[opp[i]];
			current_node.uy += ey[opp[i]]*current_node.f[opp[i]];
		}
	
		current_node.ux = current_node.ux/current_node.rho;
		current_node.uy = current_node.uy/current_node.rho;
	
		// APPLY BOUNDARY CONDITION
		if (boundary_type>0) boundary_conditions[boundary_type-1](&current_node, &boundary_value);
	
		// COLLISION
		collision_functions[collision_type](&current_node, opp, ex, ey, omega, tau, B);

		// COALESCED WRITE
		__syncthreads();
		for(int i=0;i<Q;i++)
		{
			i2d = (x + y*length.x)+i*(domain_size);
			lattice->f_curr[i2d] = current_node.f[i];
		}

		// STORE MACROS IF REQUIRED
		if (store_macros)
		{
				i2d = (x + y*length.x);
				lattice->ux[i2d] = current_node.ux;
				lattice->uy[i2d] = current_node.uy;
				lattice->u[i2d] = sqrt(current_node.ux*current_node.ux+current_node.uy*current_node.uy);
				lattice->rho[i2d] = current_node.rho;
		} 
	}
}

#endif
