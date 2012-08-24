#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include "macros.cu"

// Define a struct which represents the D2Q9 lattice
typedef struct 
{
// Population distributions
	double **f;
} Lattice;

typedef struct
{
	int *macro_bc;
	int *micro_bc;
	double *geometry;
	double **force;
	double **u;
	double *rho;
} Domain;

typedef struct
{
	int opp[Q];
	int e[DIM][Q];
	double omega[Q];
	double tau;
	double tau_mrt[Q];
	int h;
	int dt;
	int length[DIM];
	bool forcing;
	bool macro_bc;
	bool micro_bc;
	int collision_type;
	int init_type;
	double residual[NUM_RESIDS];
	double tolerance;
	double c_smag;
	double M[Q][Q];
	double M_inv[Q][Q];
} DomainConstant;

typedef struct
{
	int ixd;
	double f[Q];
	double rho;
	double u[DIM];
	double F[DIM];
	double B;
	double c_smag;
	int coord[DIM];
} Node;

typedef struct
{
	int max;
	int plot;
	int screen;
	int steady_check;
} Timing;

typedef struct
{
	bool u[DIM];
	bool rho;
	bool pressure;
	bool interactive;
	int screen_node[DIM];
} OutputController;

typedef struct
{
	char name[STR_LENGTH];
	char domain_fname[STR_LENGTH];
	char output_fname[STR_LENGTH];
} ProjectStrings;

// Solver function pointers for boundary conditions and collisions
typedef void (*micro_condition) (Node *, Lattice *);
typedef void (*macro_condition) (Node *, Domain *);
typedef void (*collision) (Node *, double *);

#endif