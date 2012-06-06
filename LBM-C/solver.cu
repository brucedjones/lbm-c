#ifndef SOLVER
#define SOLVER

#include "solver.cuh"
#include "collision.cu"

#ifdef D2Q9
	#include "boundary_conditions/d2q9_boundary.cu"
#endif

#ifdef D3Q15
	#include "boundary_conditions/d3q15_boundary.cu"
#endif



__global__ void iterate_kernel (Lattice *lattice, DomainArray *domain_arrays, bool store_macros)
{
	// Declare Variables
	double B, boundary_value;
	int ixd, target_ixd, coord[DIM], domain_size, boundary_type;
	int i, d;
	Node current_node;

	// Initialise variables
	current_node.rho = 0; current_node.u[0] = 0; current_node.u[1] = 0;
	#if DIM > 2
		current_node.u[2] = 0;
	#endif
	
	// Compute coordinates
	coord[0] = (blockDim.x*blockIdx.x)+threadIdx.x;
	coord[1] = (blockDim.y*blockIdx.y)+threadIdx.y;
	#if DIM>2
		coord[2] = (blockDim.z*blockIdx.z)+threadIdx.z;
	#endif
		
	// Load domain configuration
	double tau = domain_constants.tau;

	#if DIM > 2
		ixd = (coord[0] + coord[1]*domain_constants.length[0] + coord[2]*domain_constants.length[0]*domain_constants.length[1]);
		domain_size = domain_constants.length[0]*domain_constants.length[1]*domain_constants.length[2];
	#else
		ixd = (coord[0] + coord[1]*domain_constants.length[0]);
		domain_size = domain_constants.length[0]*domain_constants.length[1];
	#endif
	current_node.c_smag = domain_constants.c_smag;
	
	// Out-of-bounds check
	#if DIM > 2
		if(coord[0]<domain_constants.length[0] && coord[1]<domain_constants.length[1] && coord[2]<domain_constants.length[2])
	#else
		if(coord[0]<domain_constants.length[0] && coord[1]<domain_constants.length[1])
	#endif
	{
		// Set collision type and optional forces
		// The type specified in domain_constants must be multiplied by two to match the listing
		// order in the collision_functions array, an additional 1 is added to the collision type
		// to specify a collision with guo body forces
		int collision_modifier = 0;
		if(domain_constants.forcing==true)
		{
			#pragma unroll
			for (d=0;d<DIM;d++)
			{
				current_node.F[d] = domain_arrays->force[d][ixd];
				if(current_node.F[d]>0) collision_modifier = 1;
			}
		}

		int collision_type = (domain_constants.collision_type*2)+collision_modifier;

		// Load boundary condition
		if(domain_constants.zhou_he==true)
		{
			boundary_type = domain_arrays->boundary_type[ixd];
			boundary_value = domain_arrays->boundary_value[ixd];
		}	else {
			boundary_type = 0;
		}
	
		// Load Geometry
		B = domain_arrays->geometry[ixd];
		if(B==1) collision_type = 4;
	
		// STREAMING - UNCOALESCED READ
		int target_coord[DIM];
		#pragma unroll
		for(i = 0; i<Q; i++)
		{
			#pragma unroll
			for(d=0; d<DIM; d++)
			{
				target_coord[d] = coord[d]+domain_constants.e[d][i];
				if(target_coord[d]>(domain_constants.length[d]-1)) target_coord[d] = 0; if(target_coord[d]<0) target_coord[d] = domain_constants.length[d]-1;
			}

			#if DIM > 2
				target_ixd = (target_coord[0] + target_coord[1]*domain_constants.length[0] + target_coord[2]*domain_constants.length[0]*domain_constants.length[1]);
			#else
				target_ixd = (target_coord[0] + target_coord[1]*domain_constants.length[0]);
			#endif
				
			
			// UNCOALESCED READ
			current_node.f[domain_constants.opp[i]] = lattice->f_prev[domain_constants.opp[i]][target_ixd];
			current_node.rho += current_node.f[domain_constants.opp[i]];
			#pragma unroll
			for (d = 0; d<DIM; d++)
			{
				current_node.u[d] += domain_constants.e[d][domain_constants.opp[i]]*current_node.f[domain_constants.opp[i]];
			}
		}
		
		#pragma unroll
		for (d = 0; d<DIM; d++)
		{
			current_node.u[d] = current_node.u[d]/current_node.rho;
		}
		// APPLY BOUNDARY CONDITION
		if (boundary_type>0) boundary_conditions[boundary_type-1](&current_node, &boundary_value);	
		// COLLISION
		collision_functions[collision_type](&current_node, &tau, &B);
		
		// COALESCED WRITE
		__syncthreads();
		#pragma unroll
		for(int i=0;i<Q;i++)
		{
			lattice->f_curr[i][ixd] = current_node.f[i];
		}

		// STORE MACROS IF REQUIRED
		if (store_macros)
		{
			#pragma unroll
			for (d = 0; d<DIM; d++)
			{
				lattice->u[d][ixd] = current_node.u[d];
			}
			lattice->rho[ixd] = current_node.rho;
		} 
	}
}

#endif
