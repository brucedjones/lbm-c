#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include "macros.cu"

// Define a struct which represents the D2Q9 lattice
typedef struct 
{
// Population distributions
	double *f_prev;
	double *f_curr;
	double *ux;
	double *uy;
	double *u;
	double *rho;
} Lattice;

typedef struct
{
	double tau;
	int2 length;
	int b_o[5];
	int *boundary_type;
	double *boundary_value;
	double *geometry;
} Domain;

typedef struct
{
	double f[Q];
	double rho;
	double ux;
	double uy;
} Node;

// Boundary condition function pointers
typedef Node (*boundary_condition) (Node, double);

#endif