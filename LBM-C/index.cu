#ifndef KERNEL
#define KERNEL
////////////////////////////////////////////////////////////////////////////////
//
// LBM-C
// A lattice Boltzmann fluid flow solver written using CUDA
//
// Copyright (C) 2011  Bruce Jones
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
////////////////////////////////////////////////////////////////////////////////
//
// D2Q9 Lattice configuration:
//
//       6   2   5
//        \  |  /
//         \ | /
//          \|/
//       3---0---1
//          /|\
//         / | \
//        /  |  \
//       7   4   8
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "data_types.cuh"
#include "macros.cu"
#include "solver.cu"
#include "index.cuh"

// Include THRUST libraries
#include <thrust/device_vector.h>
#include <thrust/transform_reduce.h>

// DEVICE VARIABLE DECLARATION
Lattice *lattice_device;
Domain *domain_device;
double *f_1_device, *f_2_device, *rho_device, *ux_device, *uy_device, *u_device, *boundary_value_device, *geometry_device; 
int *boundary_type_device;

// HOST VARIABLE DECLARATION
Lattice *lattice_host, *lattice_device_prototype;
Domain *domain_host;
double *f_host, *rho_host, *ux_host, *uy_host, *u_host, *boundary_value_host, *geometry_host;
int *boundary_type_host;

// SCALAR DECLARATION (PLATFORM AGNOSTIC)
double tau, residual;
double tolerance;
int domain_size, l_b_o, maxT, saveT, steadyT;
int3 length;
bool store_macros = false;

int main(int argc, char **argv)
{
	//tolerance = 0.00000001;

	// Get available memory on graphics card before allocation
	size_t freeMemory_before;
	size_t totalMemory_before;
	cudaMemGetInfo(&freeMemory_before, &totalMemory_before);
	
	// Initialise memory for LBM model
	setup();

	// Get available memory on graphics card after allocation
	size_t freeMemory_after;
	size_t totalMemory_after;
	cudaMemGetInfo(&freeMemory_after, &totalMemory_after);

	// Report program memory usage
	printf("Total Device Memory:	%luMb\n", (unsigned long) totalMemory_after / 1024 / 1024);
	printf("Total Availabe Memory:	%luMb\n", (unsigned long) freeMemory_before / 1024 / 1024);
	printf("Memory Used:		%luMb\n\n", (unsigned long) (freeMemory_before-freeMemory_after) / 1024 / 1024);

	// Report domain configuration
	printf("Length.x:		%d\n", domain_host->length.x);
	printf("Length.y:		%d\n", domain_host->length.y);
	printf("Relaxation Time (Tau):	%f\n", domain_host->tau);
	printf("\nPress return to continue...");
	getchar();

	residual = 0;
	output_macros(-1);

	// Get current clock cycle number
	clock_t t1=clock();

	for(int i = 0; i<maxT; i++)
	{
		//iterate();
		if(i%saveT == 0 && steadyT>0 && i%steadyT)
		{
			store_macros = true;
			iterate();
			output_macros(i);
			residual = error_RMS(u_device,domain_size);
			if(residual<tolerance) break;
			store_macros = false;
		} else if (i%saveT==0)
		{
			store_macros = true;
			iterate();
			output_macros(i);
			store_macros = false;
		} else if(steadyT>0 && i%steadyT)
		{
			store_macros = true;
			iterate();
			cudasafe(cudaMemcpy(u_host, u_device, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");
			residual = error_RMS(u_device,domain_size);
			if(residual<tolerance) break;
			store_macros = false;
		} else{
			iterate();
		}
	}

	// Get current clock cycle number
	clock_t t2=clock();
	// Compare and report global execution time
	double cputime = ((double)t2-(double)t1)/(double)CLOCKS_PER_SEC;
	printf("\n\nTotal Run Time: %fs",cputime);
	printf("\nPress return to finish");
	getchar();


}

// ALLOCATES MEMORY ON THE HOST
void allocate_memory_host(void)
{
	// ALLOCATE ARRAY AND STRUCT MEMORY ON HOST
	// STRUCTS:
	lattice_host = (Lattice *)malloc(sizeof(Lattice));
	domain_host = (Domain *)malloc(sizeof(Domain));
	// ARRAYS:
	boundary_type_host = (int *)malloc(domain_size*sizeof(int));
	boundary_value_host = (double *)malloc(domain_size*sizeof(double));
	geometry_host = (double *)malloc(domain_size*sizeof(double));
	f_host = (double *)malloc(domain_size*Q*sizeof(double));
	rho_host = (double *)malloc(domain_size*sizeof(double));
	ux_host = (double *)malloc(domain_size*sizeof(double));
	uy_host = (double *)malloc(domain_size*sizeof(double));
	u_host = (double *)malloc(domain_size*sizeof(double));
}

// ALLOCATES MEMORY ON THE DEVICE
void allocate_memory_device(void)
{
	// ALLOCATE ARRAY AND STRUCT MEMORY ON DEVICE
	// STRUCTS:
	cudasafe(cudaMalloc((void **)&lattice_device,sizeof(Lattice)), "Allocate Memory: lattice_device");
	cudasafe(cudaMalloc((void **)&domain_device,sizeof(Domain)), "Allocate Memory: control_device");
	// ARRAYS:
	cudasafe(cudaMalloc((void **)&f_1_device,domain_size*Q*sizeof(double)), "Allocate Memory: f_1_device");
	cudasafe(cudaMalloc((void **)&f_2_device,domain_size*Q*sizeof(double)), "Allocate Memory: f_2_device");
	cudasafe(cudaMalloc((void **)&rho_device,domain_size*Q*sizeof(double)), "Allocate Memory: rho_device");
	cudasafe(cudaMalloc((void **)&ux_device,domain_size*Q*sizeof(double)), "Allocate Memory: ux_device");
	cudasafe(cudaMalloc((void **)&uy_device,domain_size*Q*sizeof(double)), "Allocate Memory: uy_device");
	cudasafe(cudaMalloc((void **)&u_device,domain_size*Q*sizeof(double)), "Allocate Memory: u_device");
	cudasafe(cudaMalloc((void **)&boundary_type_device,domain_size*sizeof(int)), "Allocate Memory: boundary_type_device");
	cudasafe(cudaMalloc((void **)&boundary_value_device,domain_size*sizeof(double)), "Allocate Memory: boundary_value_device");
	cudasafe(cudaMalloc((void **)&geometry_device,domain_size*sizeof(double)), "Allocate Memory: geometry_device");

}

// READS INPUT DATA FROM FILE AND ASSEMBLES DATA INTO RELEVANT STRUCTS
void load_and_assemble_data(void)
{
	// ASSEMBLE STRUCT ON HOST: Lattice
	lattice_host->f_curr = f_host;
	lattice_host->rho = rho_host;
	lattice_host->ux = ux_host;
	lattice_host->uy = uy_host;
	lattice_host->u = u_host;

	// ASSEMBLE AND LOAD STRUCT ON HOST: Control
	// ASSEMBLE
	domain_host->boundary_type = boundary_type_host;
	domain_host->boundary_value = boundary_value_host;
	domain_host->geometry = geometry_host;
	// LOAD
	domain_host->tau = tau;
	domain_host->length.x = length.x;
	domain_host->length.y = length.y;
	
	// Boundary nodes are treated as chains of face nodes, vertex nodes and corner nodes,
	// the length of each of these chains is a function of domain dimensions and is calculated
	// here.
	domain_host->b_o[0] = (length.y-2); // X-
	domain_host->b_o[1] = domain_host->b_o[0]+(length.y-2); // X+
	domain_host->b_o[2] = domain_host->b_o[1]+(length.x-2); // Y-
	domain_host->b_o[3] = domain_host->b_o[2]+(length.x-2); // Y+
	domain_host->b_o[4] = domain_host->b_o[3]+4;
	l_b_o = domain_host->b_o[4];

	// ASSEMBLE STRUCT ON DEVICE: Lattice
	lattice_device_prototype = (Lattice *)malloc(sizeof(Lattice));
	lattice_device_prototype->f_curr = f_1_device;
	lattice_device_prototype->f_prev = f_2_device;
	lattice_device_prototype->rho = rho_device;
	lattice_device_prototype->ux = ux_device;
	lattice_device_prototype->uy = uy_device;
	lattice_device_prototype->u = u_device;
	cudasafe(cudaMemcpy(lattice_device, lattice_device_prototype, sizeof(Lattice),cudaMemcpyHostToDevice),"Copy Data: lattice_device");

	// ASSEMBLE AND LOAD STRUCT ON DEVICE: Control
	Domain *domain_tmp = (Domain *)malloc(sizeof(Domain));
	domain_tmp->tau = tau;
	domain_tmp->length.x = length.x;
	domain_tmp->length.y = length.y;
	domain_tmp->boundary_type = boundary_type_device;
	domain_tmp->boundary_value = boundary_value_device;
	domain_tmp->geometry = geometry_device;
	cudasafe(cudaMemcpy(domain_device, domain_tmp, sizeof(Domain),cudaMemcpyHostToDevice),"Copy Data: control_device");
	cudasafe(cudaMemcpy(&domain_device->b_o, &domain_host->b_o, sizeof(int)*5,cudaMemcpyHostToDevice),"Copy Data: b_o");
}

// CALCULATES AND LOADS A CONSTANT DENSITY ZERO VELOCITY INITIAL CONDITION FOR THE DOMAIN
void load_static_IC(void)
{
	int index_i;
	double omega[Q];
	LOAD_OMEGA(omega);
	for(int i=0;i<Q;i++)
	{
		for(int index=0;index<(domain_size);index++)
		{
			index_i = index+i*(domain_size);
			lattice_host->f_curr[index_i] = 1.f*omega[i];
		}
	}
	cudasafe(cudaMemcpy(f_2_device, f_host, sizeof(double)*Q*domain_size,cudaMemcpyHostToDevice),"Copy Data: Initial Condition");
}

// EXECUTES ALL ROUTINES REQUIRED FOR THE MODEL SET UP
void setup(void)
{
	// Set cuda device to use
	cudaSetDevice(0);

	// Read domain configuration
	FILE * input_file;
    input_file = fopen ("input.dat","r");
	int IC_type, i2d;
	//IC_type = 0;
	fscanf(input_file,"%d %d %lf %d %d %d %lf %d\n", &length.x, &length.y, &tau, &saveT, &maxT, &steadyT, &tolerance, &IC_type);
	domain_size = length.x*length.y;
	allocate_memory_host();
	allocate_memory_device();
	load_and_assemble_data();
	if (IC_type == 0) load_static_IC();
	for(int j = 0; j<length.y; j++)
	{
		for(int i = 0; i<length.x; i++)
		{
			i2d = i + j*length.x;
			fscanf(input_file,"%d %lf\n", &domain_host->boundary_type[i2d], &domain_host->boundary_value[i2d]);
		}
	}

	for(int j = 0; j<length.y; j++)
	{
		for(int i = 0; i<length.x; i++)
		{
			i2d = i + j*length.x;
			fscanf(input_file,"%f\n", &domain_host->geometry[i2d]);
		}
	}

	cudasafe(cudaMemcpy(boundary_type_device, boundary_type_host, sizeof(int)*domain_size,cudaMemcpyHostToDevice),"Copy Data: boundary_type_device");
	cudasafe(cudaMemcpy(boundary_value_device, boundary_value_host, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Copy Data: boundary_value_device");
	cudasafe(cudaMemcpy(geometry_device, geometry_host, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Copy Data: geometry_device");

}

// ERROR CHECKING FOR MEMORY ALLOCATION
void cudasafe( cudaError_t error, char* message)
{
   if(error!=cudaSuccess) { fprintf(stderr,"ERROR: %s : %i\n",message,error); exit(-1); }
}

// ERROR CHECKING FOR KERNEL EXECUTION
void Check_CUDA_Error(const char *message)
{
   cudaError_t error = cudaGetLastError();
   if(error!=cudaSuccess) {
      fprintf(stderr,"ERROR: %s: %s\n", message, cudaGetErrorString(error) );
      exit(-1);
   }                         
}

// COPIES f_i DATA FROM DEVICE TO HOST AND COMPUTERS MACROSCOPIC VALUES ON HOST, THIS DATA
// IS THEN WRITTEN TO THE OUTPUT FILE
//
// Note:	A computationally more efficient implementation would compute macroscopic
//			value's on the gpu and then just copy that data, this would however consume
//			more memory
void output_macros(int time)
{
	// Copy data from device to host
	cudasafe(cudaMemcpy(rho_host, rho_device, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - rho");
	cudasafe(cudaMemcpy(ux_host, ux_device, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - ux");
	cudasafe(cudaMemcpy(uy_host, uy_device, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - uy");
	cudasafe(cudaMemcpy(u_host, u_device, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");
	
// Assemble formatted filename	
	char fname[19];
	FILE *file;
	sprintf(fname, "results_%i.dat", time);
// Open file
	file = fopen(fname,"w");
// Write File Header	
	fprintf(file,"TITLE=\"2D Poiseuille Flow\"\nVARIABLES= \"X\", \"Y\", \"rho\", \"uX\", \"uY\", \"u\"");//\nDATASETAUXDATA ComputerTime=\"%lus\"\nDATASETAUXDATA DeviceMemoryUsed=\"%luMb\"",cputime, mem);
// Write Zone Header
	// note: nx and ny values are not in the "correct" order in the zone header, errors occur when loading the data in tecplot
	// if the "correct" order is used
	fprintf(file,"\nZONE T=\"2D Poiseuille Flow at time = %i\", I=%i, J=%i, DATAPACKING=POINT, SOLUTIONTIME=%i", time,length.x,length.y,time);
// Loop over all nodes to calculate and print nodal macroscopic values to file, output some feedback data to console
	int i2d;
	for (int y=0;y<length.y;y++){
		for (int x=0;x<length.x;x++){
			
			i2d = x+y*length.x;

			// Impose zero velocity on bounceback nodes
			if(domain_host->geometry[i2d] == 1)
			{
				lattice_host->ux[i2d] = 0;
				lattice_host->uy[i2d] = 0;
				lattice_host->u[i2d] = 0;
			}
			// Write to files
			fprintf(file,"\n%i %i %lf %lf %lf %lf", x, y, lattice_host->rho[i2d], lattice_host->ux[i2d], lattice_host->uy[i2d], lattice_host->u[i2d]);
			// Output reference information to console
			if (y==length.y/2 && x == 0) {printf("\n time = %i; rho = %lf; uX = %lf; uY = %lf, resid = %g", time, lattice_host->rho[i2d], lattice_host->ux[i2d], lattice_host->uy[i2d], residual);}
			//if (y==length.y/2 && x == 0) {printf("\n time = %i; resid = %e", time, residual);}
		}
	}
	// Close file
	fclose(file);
}

// CONFIGURES THE KERNEL CONFIGURATION AND LAUNCHES KERNEL A KERNEL BOTH FOR THE BOUNDARY NODES
// AND THE BULK NODES
//
// Note:	The "all" kernel operates cleanly on both bulk and boundary nodes and may be used
//			instead, though its use is inefficient.
void iterate(void)
{
	// GRID AND BLOCK DEFINITIONS CAN BE CALCULATED BEFORE ITERATE
	// DEFINE BULK GRID AND BLOCK
	dim3 Db_bulk = dim3(length.x-2,1,1);
    dim3 Dg_bulk = dim3(length.y-2,1,1);
	// DEFINE BOUNDARY GRID AND BLOCK
	int boundary_amount = l_b_o;
	int boundary_grid=(int)(boundary_amount/BLOCK_SIZE);
	int boundary_leftover=(boundary_amount%BLOCK_SIZE);

	// ITERATE ONCE
	iterate_bulk_kernel<<<Dg_bulk, Db_bulk>>>(lattice_device, domain_device, store_macros);
	Check_CUDA_Error("Kernel \"iterate_bulk 1\" Execution Failed!");  
	iterate_boundary_kernel<<<boundary_grid,BLOCK_SIZE>>>(lattice_device,domain_device,0, store_macros);
	Check_CUDA_Error("Kernel \"iterate_boundary 1a\" Execution Failed!");  
	if(boundary_leftover)
		iterate_boundary_kernel<<<1,boundary_leftover>>>(lattice_device,domain_device,boundary_amount-boundary_leftover, store_macros);
	Check_CUDA_Error("Kernel \"iterate_boundary 1b\" Execution Failed!");

	// SWAP CURR AND PREV LATTICE POINTERS READY FOR NEXT ITER
	cudasafe(cudaMemcpy(lattice_device_prototype, lattice_device, sizeof(Lattice),cudaMemcpyDeviceToHost),"Copy Data: Device Lattice Pointers");
	double *tmp_1 = lattice_device_prototype->f_prev;
	double *tmp_2 = lattice_device_prototype->f_curr;
	lattice_device_prototype->f_curr = tmp_1;
	lattice_device_prototype->f_prev = tmp_2;
	cudasafe(cudaMemcpy(lattice_device, lattice_device_prototype, sizeof(Lattice),cudaMemcpyHostToDevice),"Copy Data: Device Lattice Pointers");
}
/*
void iterate(void)
{
	// GRID AND BLOCK DEFINITIONS CAN BE CALCULATED BEFORE ITERATE
	// DEFINE all GRID AND BLOCK
	dim3 Db_all = dim3(length.x-2,1,1);
    dim3 Dg_all = dim3(length.y-2,length.z-2,1);
	// DEFINE all GRID AND BLOCK
	int all_amount = l_b_o;
	int all_grid=(int)(all_amount/BLOCK_SIZE);
	int all_leftover=(all_amount%BLOCK_SIZE);

	// ITERATE ONCE
	iterate_all_kernel<<<Dg_all, Db_all>>>(lattice_1_device,lattice_2_device,domain_device,0,1);
	Check_CUDA_Error("Kernel \"iterate_all 1\" Execution Failed!");  
	iterate_all_kernel<<<all_grid,BLOCK_SIZE>>>(lattice_1_device,lattice_2_device,domain_device,0,0);
	Check_CUDA_Error("Kernel \"iterate_all 1a\" Execution Failed!");  
	if(all_leftover)
		iterate_all_kernel<<<1,all_leftover>>>(lattice_1_device,lattice_2_device,domain_device,all_amount-all_leftover,0);
	Check_CUDA_Error("Kernel \"iterate_all 1b\" Execution Failed!");

	// SWAP LATTICES AND ITERATE AGAIN
	iterate_all_kernel<<<Dg_all, Db_all>>>(lattice_2_device,lattice_1_device,domain_device,0,1);
	Check_CUDA_Error("Kernel \"iterate_all 2\" Execution Failed!");  
	iterate_all_kernel<<<all_grid,BLOCK_SIZE>>>(lattice_2_device,lattice_1_device,domain_device,0,0);
	Check_CUDA_Error("Kernel \"iterate_all 1a\" Execution Failed!");  
	if(all_leftover)
		iterate_all_kernel<<<1,all_leftover>>>(lattice_2_device,lattice_1_device,domain_device,all_amount-all_leftover,0);
	Check_CUDA_Error("Kernel \"iterate_all 1b\" Execution Failed!");  
}*/

// square<T> computes the square of a number f(x) -> x*x

template <typename T>
struct square
{
    __host__ __device__
        T operator()(const T& x) const { 
            return x * x;
        }
};

double current_RMS(double *device_var, int var_size)
{
	// wrap raw pointer with a device_ptr for thrust compatibility
	thrust::device_ptr<double> dev_ptr(device_var);

	// setup arguments for thrust transformation to square array elements then execute plus reduction
    square<double>        unary_op;
    thrust::plus<double> binary_op;
    double init = 0;

	// Compute RMS value
	double sum = thrust::transform_reduce(dev_ptr, dev_ptr+var_size, unary_op, init, binary_op);

	double curr_RMS = sqrt(sum/var_size);

	return curr_RMS;
}

double prev_RMS = 0;

double error_RMS(double *device_var, int var_size)
{
	double curr_RMS = current_RMS(device_var, var_size);
	double tmp = abs(curr_RMS-prev_RMS);

	prev_RMS = curr_RMS;

	return tmp;
}

#endif