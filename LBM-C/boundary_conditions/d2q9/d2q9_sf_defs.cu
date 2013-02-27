#ifndef D2Q9_SF_DEFS
#define D2Q9_SF_DEFS

__device__ __noinline__ void sf_x(Node *current_node, Lattice *lattice)
{
	// Find target index
	int target_ixd = (current_node->coord[0]+1)+current_node->coord[1]*domain_constants.length[0];

	// unknowns: 1,5,8
	current_node->f[1] = lattice->f[1][target_ixd];
	current_node->f[5] = lattice->f[5][target_ixd];
	current_node->f[8] = lattice->f[8][target_ixd];
}

__device__ __noinline__ void sf_X(Node *current_node, Lattice *lattice)
{
	// Find target index
	int target_ixd = (current_node->coord[0]-1)+current_node->coord[1]*domain_constants.length[0];

	// unknowns: 3, 6, 7
	current_node->f[3] = lattice->f[3][target_ixd];
	current_node->f[6] = lattice->f[6][target_ixd];
	current_node->f[7] = lattice->f[7][target_ixd];
}

__device__ __noinline__ void sf_y(Node *current_node, Lattice *lattice)
{
	// Find target index
	int target_ixd = current_node->coord[0]+(current_node->coord[1]+1)*domain_constants.length[0];

	// unknowns: 2, 5, 6
	current_node->f[2] = lattice->f[2][target_ixd];
	current_node->f[5] = lattice->f[5][target_ixd];
	current_node->f[6] = lattice->f[6][target_ixd];
}

__device__ __noinline__ void sf_Y(Node *current_node, Lattice *lattice)
{
	// Find target index
	int target_ixd = current_node->coord[0]+(current_node->coord[1]-1)*domain_constants.length[0];

	// unknowns: 4, 7, 8
	current_node->f[4] = lattice->f[4][target_ixd];
	current_node->f[7] = lattice->f[7][target_ixd];
	current_node->f[8] = lattice->f[8][target_ixd];
}

#endif
