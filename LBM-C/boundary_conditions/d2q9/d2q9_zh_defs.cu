#ifndef D2Q9_ZH_DEFS
#define D2Q9_ZH_DEFS

__device__ __noinline__ void zh_pressure_x(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = 1.0-((1.0/current_node->rho)*(current_node->f[0]+current_node->f[2]+current_node->f[4]+2.0*(current_node->f[3]+current_node->f[6]+current_node->f[7])));
	current_node->u[1] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[1] = current_node->f[3] + ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[5] = current_node->f[7] - ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) + ((1.0/6.0)*current_node->rho*current_node->u[0]);
	current_node->f[8] = current_node->f[6] + ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) + ((1.0/6.0)*current_node->rho*current_node->u[0]);

}

__device__ __noinline__ void zh_pressure_X(Node *current_node, Lattice *lattice)
{
	// COMPUTE MACROS
	current_node->u[0] = -1.0+((1.0/current_node->rho)*(current_node->f[0]+current_node->f[2]+current_node->f[4]+2.0*(current_node->f[1]+current_node->f[5]+current_node->f[8])));
	current_node->u[1] = 0.0;

	// COMPUTE UNKNOWN f's
	current_node->f[3] = current_node->f[1] - ((2.0/3.0)*current_node->rho*current_node->u[0]);
	current_node->f[6] = current_node->f[8] - ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) - ((1.0/6.0)*current_node->rho*current_node->u[0]);
	current_node->f[7] = current_node->f[5] + ((1.0/2.0)*(current_node->f[2]-current_node->f[4])) - ((1.0/6.0)*current_node->rho*current_node->u[0]);
}

#endif
