#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include "macros.cu"

// Define a struct which represents the D2Q9 lattice
typedef struct 
{
// Population distributions
	float *f;
} Lattice;

// Define a struct which wraps arrays for output data
typedef struct 
{
	float *rho;
	float *ux;
	float *uy;
	float *u;
} Output;

typedef struct
{
	float tau;
	int2 length;
	int b_o[5];
	int *boundary_type;
	float *boundary_value;
	float *geometry
} Domain;

typedef struct
{
	float f[Q];
	float rho;
	float ux;
	float uy;
} Node;

typedef Node (*boundary_condition) (Node, float);

#endif