#ifndef D2Q9_BOUNDARY
#define D2Q9_BOUNDARY

#include "data_types.cuh"
#include "d2q9_boundary.cuh"

__device__ inline Node zh_pressure_x(Node input, float rho_boundary)
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

__device__ inline Node zh_pressure_X(Node input, float rho_boundary)
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

__device__ inline Node zh_pressure_edge(Node input, float rho_boundary, int vector_order[8], int direction)
{
	// Applies Zhou/He pressure boundary condition for edge nodes
	// vector_order:
	//				  /	1 - Normal to boundary facing inwards
	// Unkown		-|	2 - adjacent to boundary anticlockwise
	//				  \	3 - adjacent to boundary clockwise
	//						
	//				  /	4 - Perpendicular to normal anticlockwise
	//				 |	5 - Perpendicular to normal clockwise
	// Known		-|	6 - opposite to 1
	//				 |	7 - opposite to 2
	//				  \	8 - opposite to 3
	
	Node output; //= input;

	// COPY KNOWN f's
	output.f[0] = input.f[0];
	output.f[vector_order[4]] = input.f[vector_order[4]];
	output.f[vector_order[5]] = input.f[vector_order[5]];
	output.f[vector_order[6]] = input.f[vector_order[6]];
	output.f[vector_order[7]] = input.f[vector_order[7]];
	output.f[vector_order[8]] = input.f[vector_order[8]];

	// COMPUTE MACROS
	output.rho = rho_boundary;
	output.ux = direction*1.f-direction*((1.f/output.rho)*(input.f[0]+input.f[vector_order[4]]+input.f[vector_order[5]]+2.f*(input.f[vector_order[6]]+input.f[vector_order[7]]+input.f[vector_order[8]])));
	output.uy = 0.f;

	// COMPUTE UNKNOWN f's
	output.f[vector_order[1]] = input.f[vector_order[6]] + direction*((2.f/3.f)*output.rho*output.ux);
	output.f[vector_order[2]] = input.f[vector_order[7]] + direction*((1.f/2.f)*(input.f[vector_order[4]]-input.f[vector_order[5]])) + direction*((1.f/6.f)*output.rho*output.ux);
	output.f[vector_order[3]] = input.f[vector_order[8]] - direction*((1.f/2.f)*(input.f[vector_order[4]]-input.f[vector_order[5]])) + direction*((1.f/6.f)*output.rho*output.ux);
	
	return output;
}

#endif