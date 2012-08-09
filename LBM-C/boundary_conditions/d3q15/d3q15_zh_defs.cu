#ifndef D3Q15_ZH_DEFS
#define D3Q15_ZH_DEFS

__device__ __noinline__ void zh_pressure_x(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 1.0-((1.0/current_node->rho)*(current_node->f[0]+current_node->f[3]+current_node->f[4]+current_node->f[5]+current_node->f[6]+2.0*(current_node->f[2]+current_node->f[14]+current_node->f[10]+current_node->f[12]+current_node->f[8])));
	current_node->u[1] = 0.0;
	current_node->u[2] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[1] = current_node->f[2] + ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[7] = current_node->f[14] + ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[3]-current_node->f[4]) +(current_node->f[5]-current_node->f[6])));
	current_node->f[11] = current_node->f[10] + ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[3]-current_node->f[4]) - (current_node->f[5]-current_node->f[6])));
	current_node->f[9] = current_node->f[12] + ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[3]-current_node->f[4]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[13] = current_node->f[8] + ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[3]-current_node->f[4]) - (current_node->f[5]-current_node->f[6])));
}

__device__ __noinline__ void zh_pressure_X(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = -1.0+((1.0/current_node->rho)*(current_node->f[0]+current_node->f[3]+current_node->f[4]+current_node->f[5]+current_node->f[6]+2.0*(current_node->f[1]+current_node->f[7]+current_node->f[11]+current_node->f[9]+current_node->f[13])));
	current_node->u[1] = 0.0;
	current_node->u[2] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[2] = current_node->f[1] - ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[14] = current_node->f[7] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[3]-current_node->f[4]) - (current_node->f[5]-current_node->f[6])));
	current_node->f[10] = current_node->f[11] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[3]-current_node->f[4]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[12] = current_node->f[9] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[3]-current_node->f[4]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[8] = current_node->f[13] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[3]-current_node->f[4]) - (current_node->f[5]-current_node->f[6])));
}

__device__ __noinline__ void zh_pressure_y(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 0.0;
	current_node->u[1] = 1.0-((1.0/current_node->rho)*(current_node->f[0]+current_node->f[1]+current_node->f[2]+current_node->f[5]+current_node->f[6]+2.0*(current_node->f[4]+current_node->f[14]+current_node->f[10]+current_node->f[9]+current_node->f[13])));
	current_node->u[2] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[3] = current_node->f[4] + ((2.0/3.0)*current_node->rho*current_node->u[1]);
	current_node->f[7] = current_node->f[14] + ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[11] = current_node->f[10] + ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) - (current_node->f[5]-current_node->f[6])));
	current_node->f[12] = current_node->f[9] + ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) - (current_node->f[5]-current_node->f[6])));
	current_node->f[8] = current_node->f[13] + ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) + (current_node->f[5]-current_node->f[6])));
}

__device__ __noinline__ void zh_pressure_Y(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 0.0;
	current_node->u[1] = -1.0+((1.0/current_node->rho)*(current_node->f[0]+current_node->f[1]+current_node->f[2]+current_node->f[5]+current_node->f[6]+2.0*(current_node->f[3]+current_node->f[7]+current_node->f[11]+current_node->f[12]+current_node->f[8])));
	current_node->u[2] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[4] = current_node->f[3] - ((2.0/3.0)*current_node->rho*current_node->u[1]);
	current_node->f[14] = current_node->f[7] - ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) - (current_node->f[5]-current_node->f[6])));
	current_node->f[10] = current_node->f[11] - ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[9] = current_node->f[12] - ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) + (current_node->f[5]-current_node->f[6])));
	current_node->f[13] = current_node->f[8] - ((1.0/12.0)*current_node->rho*current_node->u[1]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) - (current_node->f[5]-current_node->f[6])));
}

__device__ __noinline__ void zh_pressure_z(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 0.0;
	current_node->u[1] = 0.0;
	current_node->u[2] = 1.0-((1.0/current_node->rho)*(current_node->f[0]+current_node->f[1]+current_node->f[2]+current_node->f[3]+current_node->f[4]+2.0*(current_node->f[6]+current_node->f[14]+current_node->f[11]+current_node->f[12]+current_node->f[13])));

	// COMPUTE UNKNOWN f's
	current_node->f[5] = current_node->f[6] + ((2.0/3.0)*current_node->rho*current_node->u[2]);
	current_node->f[7] = current_node->f[14] + ((1.0/12.0)*current_node->rho*current_node->u[2]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) + (current_node->f[3]-current_node->f[4])));
	current_node->f[10] = current_node->f[11] + ((1.0/12.0)*current_node->rho*current_node->u[2]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) - (current_node->f[3]-current_node->f[4])));
	current_node->f[9] = current_node->f[12] + ((1.0/12.0)*current_node->rho*current_node->u[2]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) - (current_node->f[3]-current_node->f[4])));
	current_node->f[8] = current_node->f[13] + ((1.0/12.0)*current_node->rho*current_node->u[2]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) + (current_node->f[3]-current_node->f[4])));
}

__device__ __noinline__ void zh_pressure_Z(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 0.0;
	current_node->u[1] = 0.0;
	current_node->u[2] =  -1.0+((1.0/current_node->rho)*(current_node->f[0]+current_node->f[1]+current_node->f[2]+current_node->f[3]+current_node->f[4]+2.0*(current_node->f[5]+current_node->f[7]+current_node->f[10]+current_node->f[9]+current_node->f[8])));

	// COMPUTE UNKNOWN f's
	current_node->f[6] = current_node->f[5] - ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[14] = current_node->f[7] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) - (current_node->f[3]-current_node->f[4])));
	current_node->f[11] = current_node->f[10] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) + (current_node->f[3]-current_node->f[4])));
	current_node->f[12] = current_node->f[9] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((current_node->f[1]-current_node->f[2]) + (current_node->f[3]-current_node->f[4])));
	current_node->f[13] = current_node->f[8] - ((1.0/12.0)*current_node->rho*current_node->u[0]) - ((1.0/4.0)*((-1.0)*(current_node->f[1]-current_node->f[2]) - (current_node->f[3]-current_node->f[4])));
}
#endif
