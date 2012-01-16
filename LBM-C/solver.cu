#ifndef SOLVER
#define SOLVER

#include "solver.cuh"
#include "d2q9_boundary.cu"


__device__ boundary_condition boundary_conditions[2] = { zh_pressure_x, zh_pressure_X};

// PREFORMS ONE ITERATION OF THE LBM ON BULK NODES(NODES WHICH ARE NOT ON A DOMAIN BOUNDARY)
__global__ void iterate_bulk_kernel (Lattice *lattice, Domain *domain, bool store_macros)
{
	// Declare variables
	double f_eq,omega[Q],cu,u_sq, collision_bgk, collision_s, B;
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
	int x     = threadIdx.x+1;
    int y     = blockIdx.x+1;

	// Load domain configuration
	length.x = domain->length.x;
	length.y = domain->length.y;
	int domain_size = length.x*length.y;
	double tau = domain->tau;

	// Load geometry
	int i2d_prime = x + y*length.x;
	B = domain->geometry[i2d_prime];

	// STREAMING - Stream f's and calculate macroscopic values. Streaming occurs here as streaming represents
	// an uncoalesced memory access, uncoalesced reads are less time consuming than uncoalesced
	// writes.
	int target_x, target_y;
	for(int i = 0; i<Q; i++)
	{
		target_x = x+ex[i]; target_y = y+ey[i];

		i2d = (target_x + target_y*length.x)+opp[i]*(domain_size);
		
		// UNCOALESCED READ
		current_node.f[opp[i]] = lattice->f_prev[i2d];

		current_node.rho += current_node.f[opp[i]];
		current_node.ux += ex[opp[i]]*current_node.f[opp[i]];
		current_node.uy += ey[opp[i]]*current_node.f[opp[i]];
	}
	
	current_node.ux = current_node.ux/current_node.rho;
	current_node.uy = current_node.uy/current_node.rho;	

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
	u_sq = 1.5f*(current_node.ux*current_node.ux + current_node.uy*current_node.uy);
	for(int i=0;i<Q;i++)
	{
		i2d = (x + y*length.x)+i*(domain_size);

		cu = 3*(ex[i]*current_node.ux+ey[i]*current_node.uy);
		f_eq = current_node.rho*omega[i]*(1.f+cu+(0.5f*cu*cu)-u_sq);

		collision_bgk = (1.f/tau) * (current_node.f[i]-f_eq);
		collision_s = current_node.f[opp[i]]-current_node.f[i];

		lattice->f_curr[i2d] = current_node.f[i] - (1-B)*collision_bgk + B*collision_s;
	}
}

// PREFORMS ONE ITERATION OF THE LBM ON BOUNDARY NODES
__global__ void iterate_boundary_kernel (Lattice *lattice, Domain *domain, int offset, bool store_macros)
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
	int idx=blockIdx.x*BLOCK_SIZE+threadIdx.x+offset;
	int2 coords = compute_boundary_coords(idx, domain);
	int x = coords.x;
	int y = coords.y;

	// Load domain configuration
	length.x = domain->length.x;
	length.y = domain->length.y;
	int domain_size = length.x*length.y;
	double tau = domain->tau;
	int i2d_prime = x + y*length.x;

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
		for(int i=0;i<Q;i++)
		{
			i2d = (x + y*length.x);
			lattice->ux[i2d] = current_node.ux;
			lattice->uy[i2d] = current_node.uy;
			lattice->u[i2d] = sqrt(current_node.ux*current_node.ux+current_node.uy*current_node.uy);
			lattice->rho[i2d] = current_node.rho;
		}
	} 

	// COLLISION - COALESCED WRITE
	u_sq = 1.5f*(current_node.ux*current_node.ux + current_node.uy*current_node.uy);
	for(int i=0;i<Q;i++)
	{
		i2d = (x + y*length.x)+i*(domain_size);

		cu = 3*(ex[i]*current_node.ux+ey[i]*current_node.uy);
		f_eq = current_node.rho*omega[i]*(1.f+cu+(0.5f*cu*cu)-u_sq);

		collision_bgk = (1.f/tau) * (current_node.f[i]-f_eq);
		collision_s = current_node.f[opp[i]]-current_node.f[i];

		lattice->f_curr[i2d] = current_node.f[i] - (1-B)*collision_bgk + B*collision_s;
	}
}

// KERNEL SUPPOT FUNCTION
__device__ inline int2 compute_boundary_coords(int idx, Domain *domain)
{
	// All boundary nodes are indexed by a single index, thresholds define
	// the length of boundary edges and the way in which the x and y 
	// coords are calculated from the index

	int2 coord, length;
	int id;

	length.x = domain->length.x;
	length.y = domain->length.y;

	if(idx>=0 && idx<domain->b_o[0]) // X-
	{
		id = idx;
		coord.y = id+1;
		coord.x = 0;
	} else if (idx>=domain->b_o[0] && idx<domain->b_o[1]) // X+
	{
		id = idx-domain->b_o[0];
		coord.y = id+1;
		coord.x = length.x-1;
	} else if (idx>=domain->b_o[1] && idx<domain->b_o[2]) // Y-
	{
		id = idx-domain->b_o[1];
		coord.x = id+1;
		coord.y = 0;
	} else if (idx>=domain->b_o[2] && idx<domain->b_o[3]) // Y+
	{
		id = idx-domain->b_o[2];
		coord.x = id+1;
		coord.y = length.y-1;
	} else if (idx>=domain->b_o[3]) //CORNERS
	{
		id = idx-domain->b_o[3];
		if(id == 0)
		{
			coord.x = 0;
			coord.y = 0;
		} else if(id==1)
		{
			coord.x = 0;
			coord.y = length.y-1;
		} else if(id==2)
		{
			coord.x = length.x-1;
			coord.y = length.y-1;
		} else if(id==3)
		{
			coord.x = length.x-1;
			coord.y = 0;
		}
	}

	return coord;
}

#endif
