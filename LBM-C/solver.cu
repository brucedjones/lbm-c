#ifndef SOLVER
#define SOLVER

#include "solver.cuh"
#include "d2q9_boundary.cu"


__device__ boundary_condition boundary_conditions[2] = { zh_pressure_x, zh_pressure_X};

// PREFORMS ONE ITERATION OF THE LBM ON BOUNDARY NODES
__global__ void iterate_kernel (Lattice *lattice, Domain *domain, bool store_macros)
{
	// Declare Variables
	double f_eq, omega[Q], cu, u_sq, collision_bgk, collision_s, B;
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
	length.x = domain->length.x;
	length.y = domain->length.y;
	int domain_size = length.x*length.y;
	double tau = domain->tau;
	int i2d_prime = x + y*length.x;

	if(x<length.x && y<length.y)
	{
	
		// Load boundary condition
		int boundary_type = domain->boundary_type[i2d_prime];
		double boundary_value = domain->boundary_value[i2d_prime];
	
		// Load Geometry
		B = domain->geometry[i2d_prime];
	
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
	
		// STORE MACROS IF REQUIRED
		if (store_macros)
		{
				i2d = (x + y*length.x);
				lattice->ux[i2d] = current_node.ux;
				lattice->uy[i2d] = current_node.uy;
				lattice->u[i2d] = sqrt(current_node.ux*current_node.ux+current_node.uy*current_node.uy);
				lattice->rho[i2d] = current_node.rho;
		} 
	
		// COLLISION - COALESCED WRITE
		u_sq = 1.5*(current_node.ux*current_node.ux + current_node.uy*current_node.uy);
		for(int i=0;i<Q;i++)
		{
			i2d = (x + y*length.x)+i*(domain_size);
	
			cu = 3.0*(ex[i]*current_node.ux+ey[i]*current_node.uy);
			f_eq = current_node.rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
			collision_bgk = (1.0/tau) * (current_node.f[i]-f_eq);
			collision_s = current_node.f[opp[i]]-current_node.f[i];
	
			lattice->f_curr[i2d] = current_node.f[i] - (1-B)*collision_bgk + B*collision_s;
		}
	}
}

__global__ void iterate_forced_kernel (Lattice *lattice, Domain *domain, bool store_macros)
{
	// Declare Variables
	double f_eq, omega[Q], cu, u_sq, collision_bgk, collision_s, B;
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
	length.x = domain->length.x;
	length.y = domain->length.y;
	int domain_size = length.x*length.y;
	double tau = domain->tau;
	int i2d_prime = x + y*length.x;

	if(x<length.x && y<length.y)
	{
	
		// Load boundary condition
		int boundary_type = domain->boundary_type[i2d_prime];
		double boundary_value = domain->boundary_value[i2d_prime];
	
		// Load Geometry
		B = domain->geometry[i2d_prime];
	
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
		if (boundary_type>0) current_node = boundary_conditions[boundary_type-1](current_node, boundary_value);
	
		// STORE MACROS IF REQUIRED
		if (store_macros)
		{
				i2d = (x + y*length.x);
				lattice->ux[i2d] = current_node.ux;
				lattice->uy[i2d] = current_node.uy;
				lattice->u[i2d] = sqrt(current_node.ux*current_node.ux+current_node.uy*current_node.uy);
				lattice->rho[i2d] = current_node.rho;
		} 
	
		// COLLISION - COALESCED WRITE
		u_sq = 1.5*(current_node.ux*current_node.ux + current_node.uy*current_node.uy);
		for(int i=0;i<Q;i++)
		{
			i2d = (x + y*length.x)+i*(domain_size);
	
			cu = 3.0*(ex[i]*current_node.ux+ey[i]*current_node.uy);
			f_eq = current_node.rho*omega[i]*(1.0+cu+(0.5*cu*cu)-u_sq);
	
			collision_bgk = (1.0/tau) * (current_node.f[i]-f_eq);
			collision_s = current_node.f[opp[i]]-current_node.f[i];
	
			lattice->f_curr[i2d] = current_node.f[i] - (1-B)*collision_bgk + B*collision_s;
		}
	}
}

#endif
