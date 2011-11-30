#ifndef D2Q9_BOUNDARY
#define D2Q9_BOUNDARY

// Necessary includes
#include "macros.cu"
#include "d2q9_boundary.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "data_types.cuh"
#include "cuda_runtime.h"

__device__ __noinline__ Node zh_pressure_x(Node input, float rho_boundary)
{
	Node output; //= input;

	// COPY KNOWN f's
	output.f[0] = input.f[0];
	output.f[2] = input.f[2];
	output.f[4] = input.f[4];
	output.f[3] = input.f[3];
	output.f[6] = input.f[6];
	output.f[7] = input.f[7];

	// COMPUTE MACROS
	output.rho = rho_boundary;
	output.ux = 1.f-((1.f/output.rho)*(input.f[0]+input.f[2]+input.f[4]+2.f*(input.f[3]+input.f[6]+input.f[7])));
	output.uy = 0.f;

	// COMPUTE UNKNOWN f's
	output.f[1] = input.f[3] + ((2.f/3.f)*output.rho*output.ux);
	output.f[5] = input.f[7] - ((1.f/2.f)*(input.f[2]-input.f[4])) + ((1.f/6.f)*output.rho*output.ux);
	output.f[8] = input.f[6] + ((1.f/2.f)*(input.f[2]-input.f[4])) + ((1.f/6.f)*output.rho*output.ux);
	
	return output;
}

__device__ __noinline__ Node zh_pressure_X(Node input, float rho_boundary)
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
	output.ux = -1.f+((1.f/output.rho)*(input.f[0]+input.f[2]+input.f[4]+2.f*(input.f[1]+input.f[5]+input.f[8])));
	output.uy = 0.f;

	// COMPUTE UNKNOWN f's
	output.f[3] = input.f[1] - ((2.f/3.f)*output.rho*output.ux);
	output.f[6] = input.f[8] - ((1.f/2.f)*(input.f[2]-input.f[4])) - ((1.f/6.f)*output.rho*output.ux);
	output.f[7] = input.f[5] + ((1.f/2.f)*(input.f[2]-input.f[4])) - ((1.f/6.f)*output.rho*output.ux);
	
	return output;
}
/*
__global__ void iterate_all_kernel (Lattice *lattice_1, Lattice *lattice_2, Domain *domain, int offset, int type)
{
	int x,y,z;
	if(type==1)
	{
		x     = threadIdx.x+1;
		y     = blockIdx.x+1;
		z     = blockIdx.y+1;
	} else
	{
		int idx=blockIdx.x*BLOCK_SIZE+threadIdx.x+offset;
		int3 coords = compute_boundary_coords(idx, domain);
		x = coords.x;
		y = coords.y;
		z = coords.z;
	}

	float f_eq, f_eqb[Q],omega[Q],cu,u_sq, collision_bgk, collision_s, B;
	int i2d, ex[Q], ey[Q], ez[Q], opp[Q];
	int3 length;
	Node current_node;
	current_node.rho = 0; current_node.ux = 0; current_node.uy = 0; current_node.uz = 0;

	LOAD_EX(ex);
	LOAD_EY(ey);
	LOAD_EZ(ez);
	LOAD_OMEGA(omega);
	LOAD_OPP(opp);
	
	length.x = domain->length.x;
	length.y = domain->length.y;
	length.z = domain->length.z;
	int domain_size = length.x*length.y*length.z;

	float tau = domain->tau;

	int target_x, target_y, target_z;

	int i2d_prime = x + y*length.x + z*length.y*length.x;
	float boundary_type = floor(domain->boundary_type[i2d_prime]);
	float boundary_value = domain->boundary_value[i2d_prime];
	if (boundary_type >= 1) 
	{
		B = domain->boundary_type[i2d_prime]-boundary_type;
	} else if (boundary_type < 1) 
	{
		B = 1.f;
	} 
	//B = 0.f;


	for(int i = 0; i<Q; i++)
	{
		target_x = x+ex[i]; target_y = y+ey[i]; target_z = z+ez[i];
		//PERIODIC BOUNDARY
		if(target_x>(length.x-1)) target_x = 0; if(target_x<0) target_x = length.x-1;
		if(target_y>(length.y-1)) target_y = 0; if(target_y<0) target_y = length.y-1;
		if(target_z>(length.z-1)) target_z = 0; if(target_z<0) target_z = length.z-1;

		i2d = (target_x + target_y*length.x + target_z*length.y*length.x)+opp[i]*(domain_size);
		
		// UNCOALESCED READ
		current_node.f[opp[i]] = lattice_1->f[i2d];

		current_node.rho += current_node.f[opp[i]];
		current_node.ux += ex[opp[i]]*current_node.f[opp[i]];
		current_node.uy += ey[opp[i]]*current_node.f[opp[i]];
		current_node.uz += ez[opp[i]]*current_node.f[opp[i]];
	}
	
	current_node.ux = current_node.ux/current_node.rho;
	current_node.uy = current_node.uy/current_node.rho;
	current_node.uz = current_node.uz/current_node.rho;

	// APPLY BOUNDARY CONDITION
	if(type!=1)
	{
		if(boundary_type == 2) current_node = zh_pressure_ZY_x(current_node, boundary_value);
		if(boundary_type == 3) current_node = zh_pressure_ZY_X(current_node, boundary_value);
	}

	u_sq = 1.5f*(current_node.ux*current_node.ux + current_node.uy*current_node.uy + current_node.uz*current_node.uz);
	

	// COALESCED WRITE
	for(int i=0;i<Q;i++)
	{
		i2d = (x + y*length.x + z*length.y*length.x)+i*(domain_size);

		cu = 3*(ex[i]*current_node.ux+ey[i]*current_node.uy+ez[i]*current_node.uz);
		f_eq = current_node.rho*omega[i]*(1.f+cu+(0.5f*cu*cu)-u_sq);

		collision_bgk = (1.f/tau) * (current_node.f[i]-f_eq);
		collision_s = current_node.f[opp[i]]-current_node.f[i];

		lattice_2->f[i2d] = current_node.f[i] - (1-B)*collision_bgk + B*collision_s;
	}
}*/


#endif
