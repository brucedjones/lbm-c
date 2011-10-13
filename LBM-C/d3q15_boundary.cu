#ifndef D3Q15_BOUNDARY
#define D3Q15_BOUNDARY

#include "data_types.cuh"
#include "d3q15_boundary.cuh"

__device__ inline Node zh_pressure_ZY_x(Node input, float rho_boundary)
{
	Node output; //= input;
	output.f[0] = input.f[0];
	output.f[2] = input.f[2];
	output.f[3] = input.f[3];
	output.f[4] = input.f[4];
	output.f[5] = input.f[5];
	output.f[6] = input.f[6];
	output.f[8] = input.f[8];
	output.f[10] = input.f[10];
	output.f[12] = input.f[12];
	output.f[14] = input.f[14];

	//float rhoIn = INLET_RHO;
	// COMPUTE MACROS
	output.rho = rho_boundary;
	output.ux = 1.f-((1.f/output.rho)*(input.f[0]+input.f[3]+input.f[4]+input.f[5]+input.f[6]+2.f*(input.f[2]+input.f[8]+input.f[10]+input.f[12]+input.f[14])));
	output.uy = 0.f;
	output.uz = 0.f;
	// COMPUTE UNKNOWN f's
	output.f[1] = input.f[2] +((2.f/3.f)*output.rho*output.ux);
	output.f[7] = input.f[8] +((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*((input.f[3]-input.f[4])+(input.f[5]-input.f[6])));
	output.f[9] = input.f[10] +((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*((input.f[3]-input.f[4])-(input.f[5]-input.f[6])));
	output.f[11] = input.f[12] +((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*(-(input.f[3]-input.f[4])+(input.f[5]-input.f[6])));
	output.f[13] = input.f[14] +((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*(-(input.f[3]-input.f[4])-(input.f[5]-input.f[6])));
	
	return output;
}

__device__ inline Node zh_pressure_ZY_X(Node input, float rho_boundary)
{
	Node output; //= input;
	output.f[0] = input.f[0];
	output.f[1] = input.f[1];
	output.f[3] = input.f[3];
	output.f[4] = input.f[4];
	output.f[5] = input.f[5];
	output.f[6] = input.f[6];
	output.f[7] = input.f[7];
	output.f[9] = input.f[9];
	output.f[11] = input.f[11];
	output.f[13] = input.f[13];

	// COMPUTE MACROS
	output.rho = rho_boundary;
	output.ux = -1.f+((1.f/output.rho)*(input.f[0]+input.f[3]+input.f[4]+input.f[5]+input.f[6]+2.f*(input.f[1]+input.f[7]+input.f[9]+input.f[11]+input.f[13])));
	output.uy = 0.f;
	output.uz = 0.f;
	// COMPUTE UNKNOWN f's
	output.f[2] = input.f[1] -((2.f/3.f)*output.rho*output.ux);
	output.f[8] = input.f[7] -((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*(-(input.f[3]-input.f[4])-(input.f[5]-input.f[6])));
	output.f[10] = input.f[9] -((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*(-(input.f[3]-input.f[4])+(input.f[5]-input.f[6])));
	output.f[12] = input.f[11] -((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*((input.f[3]-input.f[4])-(input.f[5]-input.f[6])));
	output.f[14] = input.f[13] -((1.f/12.f)*output.rho*output.ux) - ((1.f/4.f)*((input.f[3]-input.f[4])+(input.f[5]-input.f[6])));
	
	return output;
}

#endif