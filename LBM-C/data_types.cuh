#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include "macros.cu"

// Define a struct which represents the D2Q9 lattice
typedef struct 
{
// Population distributions
	double **f_prev;
	double **f_curr;
	double **u;
	double *rho;
} Lattice;

typedef struct
{
	int *boundary_type;
	double *boundary_value;
	double *geometry;
	double **force;
} DomainArray;

typedef struct
{
	double tau;
	int h;
	int dt;
	int length[DIM];
	bool forcing;
	int zhou_he;
	int collision_type;
	int init_type;
	double residual;
	double tolerance;
	double c_smag
} DomainConstant;

typedef struct
{
	double f[Q];
	double rho;
	double u[DIM];
	double F[DIM];
	double B;
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
} OutputController;

typedef struct
{
	char name[STR_LENGTH];
	char domain_fname[STR_LENGTH];
	char output_fname[STR_LENGTH];
} ProjectStrings;

// Solver function pointers for boundary conditions and collisions
typedef void (*boundary_condition) (Node *, double *);
typedef void (*collision) (Node *, int *, int [DIM][Q], double *, double *, double *);

#endif