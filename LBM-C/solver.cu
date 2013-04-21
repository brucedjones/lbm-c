#ifndef SOLVER
#define SOLVER

#include "solver.cuh"
#include "collision.cu"
#include "boundary_conditions/macro_bc.cu"
#include "boundary_conditions/micro_bc.cuh"


__global__ void iterate_kernel (Lattice *lattice, Domain *domain, bool store_macros, int t)
{
	// Declare Variables
	int ixd, target_ixd, domain_size;
	int i, d, macro_bc, micro_bc;
	int ix[DIM][Q], ixd2[Q];
	Node current_node;
	
	// Compute coordinates
	current_node.coord[0] = (blockDim.x*blockIdx.x)+threadIdx.x;
	current_node.coord[1] = (blockDim.y*blockIdx.y)+threadIdx.y;
	#if DIM>2
		current_node.coord[2] = (blockDim.z*blockIdx.z)+threadIdx.z;
		ixd = (current_node.coord[0] + current_node.coord[1]*domain_constants.length[0] + current_node.coord[2]*domain_constants.length[0]*domain_constants.length[1]);
		domain_size = domain_constants.length[0]*domain_constants.length[1]*domain_constants.length[2];
	#else
		ixd = (current_node.coord[0] + current_node.coord[1]*domain_constants.length[0]);
		domain_size = domain_constants.length[0]*domain_constants.length[1];
	#endif
		current_node.ixd = ixd;
	
	// Out-of-bounds check
	#if DIM > 2
		if(current_node.coord[0]<domain_constants.length[0] && current_node.coord[1]<domain_constants.length[1] && current_node.coord[2]<domain_constants.length[2])
	#else
		if(current_node.coord[0]<domain_constants.length[0] && current_node.coord[1]<domain_constants.length[1])
	#endif
	{
		// Initialise variables
		current_node.rho = 0; current_node.u[0] = 0; current_node.u[1] = 0;
		#if DIM > 2
			current_node.u[2] = 0;
		#endif

		// Load domain configuration
		// Relaxation time:
		double tau = domain_constants.tau;
		// Smagorinsky constant:
		current_node.c_smag = domain_constants.c_smag;
		// Boundary condition:
		int macro_bc = 0;
		if(domain_constants.macro_bc==true)	macro_bc = domain->macro_bc[ixd];
		int micro_bc = 0;
		if(domain_constants.micro_bc==true)	micro_bc = domain->micro_bc[ixd];

		// Collision type:
		int collision_type = (domain_constants.collision_type*2);	// Set collision type and optional forces
		if(domain_constants.forcing==true)							// The type specified in domain_constants must be multiplied by two to match the listing
		{	
			collision_type += 1;									// order in the collision_functions array, an additional 1 is added to the collision type
			#pragma unroll											// to specify a collision with guo body forces
			for (d=0;d<DIM;d++)
			{
				current_node.F[d] = domain->force[d][ixd];
			}
		}
		// Geometry:
		current_node.B = domain->geometry[ixd];
		if(current_node.B>=1){
			collision_type = 8;
			micro_bc = 0;
			macro_bc = 0;
		}

		// COALESCED READ
		int target_coord[DIM];
		#pragma unroll
		for(i = 0; i<Q; i++)
		{
			for(d = 0; d<DIM;d++)
			{
				// Compute PDF origin ijk
				ix[d][i] = ((current_node.coord[d]-(t*domain_constants.e[d][i]))%domain_constants.length[d]);
				// Account for C implementation of modulo on negative numbers
				if (ix[d][i]<0) ix[d][i] = domain_constants.length[d]-(ix[d][i]*-1);
			}
			
			// Compute PDF origin index
			#if DIM > 2
				ixd2[i] = ix[0][i] + ix[1][i]*domain_constants.length[0] + ix[2][i]*domain_constants.length[0]*domain_constants.length[1];
			#else
				ixd2[i] = ix[0][i] + ix[1][i]*domain_constants.length[0];
			#endif
			// COALESCED READ
			current_node.f[i] = lattice->f[i][ixd2[i]];
			// CALCULATE MACROS
			current_node.rho += current_node.f[i];
			#pragma unroll
			for (d = 0; d<DIM; d++)
			{
				current_node.u[d] += domain_constants.e[d][i]*current_node.f[i];
			}
		}
		
		#pragma unroll
		for (d = 0; d<DIM; d++)
		{
			current_node.u[d] = current_node.u[d]/current_node.rho;
		}

		// Load boundary condition
		if(macro_bc>0) macro_conditions[macro_bc-1](&current_node, domain);	
		if(micro_bc>0) micro_conditions[micro_bc-1](&current_node, lattice);

		// COLLISION
		collision_functions[collision_type](&current_node, &tau);
		
		// COALESCED STREAMING WRITE
		__syncthreads();
		#pragma unroll
		for(int i=0;i<Q;i++)
		{
			lattice->f[i][ixd2[i]] = current_node.f[i];
		}

		// STORE MACROS IF REQUIRED
		if (store_macros)
		{
			#pragma unroll
			for (d = 0; d<DIM; d++)
			{
				domain->u[d][ixd] = current_node.u[d];
			}
			domain->rho[ixd] = current_node.rho;
		} 
	}
}

#endif
