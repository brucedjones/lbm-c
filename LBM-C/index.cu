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

#ifdef _WIN64
	#pragma comment(lib, "cgns/x64/lib/cgns.lib")
	#include "cgns\x64\include\cgnslib.h"
	#pragma comment(lib, "HDF5/x64/lib/hdf5.lib")
	#include "HDF5/x64/include/hdf5.h"
	#pragma comment(lib, "HDF5/x64/lib/libszip.lib")
	#include "HDF5/x64/include/szlib.h"
	#pragma comment(lib, "HDF5/x64/lib/libzlib.lib")
	#include "HDF5/x64/include/zlib.h"
#else
	#pragma comment(lib, "cgns/x86/lib/cgns.lib")
	#include "cgns\x86\include\cgnslib.h"
	#pragma comment(lib, "HDF5/x86/lib/hdf5.lib")
	#include "HDF5/x86/include/hdf5.h"
	#pragma comment(lib, "HDF5/x86/lib/libszip.lib")
	#include "HDF5/x86/include/szlib.h"
	#pragma comment(lib, "HDF5/x86/lib/libzlib.lib")
	#include "HDF5/x86/include/zlib.h"
#endif

#include <stdio.h>
#include "data_types.cuh"
#include "macros.cu"
#include "solver.cu"
#include "index.cuh"
#include "model_builder.cu"
#include "cgns/cgns_output_handler.cu"
#include "cuda_util.cu"

// Include THRUST libraries
#include <thrust/device_vector.h>
#include <thrust/transform_reduce.h>

// DEVICE VARIABLE DECLARATION
Lattice *lattice_device;
DomainArray *domain_arrays_device;
DomainConstant *domain_constants_device;
OutputController *output_controller_device;

// HOST VARIABLE DECLARATION
Lattice *lattice_host, *lattice_device_prototype;
DomainArray *domain_arrays_host;
DomainConstant *domain_constants_host;
OutputController *output_controller_host;
Timing *times;
ProjectStrings *project;
ModelBuilder model_builder;


// SCALAR DECLARATION (PLATFORM AGNOSTIC)
bool store_macros = false;

// DECLARE OUTPUT HANDLER
CGNSOutputHandler output_handler;

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
	printf("X-Length:		%d\n", domain_constants_host->length[0]);
	printf("Y-Length:		%d\n", domain_constants_host->length[1]);
	#if DIM > 2
		printf("Z-Length:		%d\n", domain_constants_host->length[2]);
	#endif
	printf("Relaxation Time (Tau):	%f\n", domain_constants_host->tau);
	printf("\nPress return to continue...");
	getchar();

	domain_constants_host->residual = 0;
	output_macros(-1);

	// Get current clock cycle number
	clock_t t1=clock();

	int domain_size=1;
	for(int d = 0; d<DIM ;d++)
	{
		domain_size = domain_size*domain_constants_host->length[d];
	}

	for(int i = 0; i<times->max; i++)
	{
		if(i%times->plot == 0 || (times->steady_check>0 && i%times->steady_check) || i%times->screen)
		{
			store_macros = true;
		}

		iterate();
		store_macros = false;

		if (i%times->plot==0) output_macros(i);
		if (i%times->steady_check==0) steady_check;

		if(i%times->plot == 0 && times->steady_check>0 && i%times->steady_check)
		{
			store_macros = true;
			iterate();
			output_macros(i);
			domain_constants_host->residual = error_RMS(lattice_device->u[0],domain_size);
			if(domain_constants_host->residual<domain_constants_host->tolerance) break;
			store_macros = false;
		} else if (i%times->plot==0)
		{
			store_macros = true;
			iterate();
			output_macros(i);
			store_macros = false;
		} else if(times->steady_check>0 && i%times->steady_check)
		{
			store_macros = true;
			iterate();
			cudasafe(cudaMemcpy(lattice_host->u[0], lattice_device->u[0], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");
			domain_constants_host->residual = error_RMS(lattice_device->u[0],domain_size);
			if(domain_constants_host->residual<domain_constants_host->tolerance) break;
			store_macros = false;
		} else
		{
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


// EXECUTES ALL ROUTINES REQUIRED FOR THE MODEL SET UP
void setup(void)
{
	// Set cuda device to use
	cudaSetDevice(0);
	cudaFuncSetCacheConfig(iterate_kernel, cudaFuncCachePreferL1);

	// Allocate container structures
	combi_malloc<Lattice>(&lattice_host, &lattice_device, sizeof(Lattice));
	combi_malloc<DomainArray>(&domain_arrays_host, &domain_arrays_device, sizeof(DomainArray));
	combi_malloc<DomainConstant>(&domain_constants_host, &domain_constants_device, sizeof(DomainConstant));
	combi_malloc<OutputController>(&output_controller_host, &output_controller_device, sizeof(OutputController));
	domain_constants_host = (DomainConstant *)malloc(sizeof(DomainConstant));
	times = (Timing *)malloc(sizeof(Timing));
	project = (ProjectStrings *)malloc(sizeof(ProjectStrings));
	lattice_device_prototype = (Lattice *)malloc(sizeof(Lattice));

	ModelBuilder tmpmb("cylinder.lbmc", lattice_host, lattice_device,
		domain_constants_host, domain_constants_device,
		domain_arrays_host, domain_arrays_device,
		output_controller_host, output_controller_device,
		times, project);
	model_builder = tmpmb;

	/*model_builder.get_model(lattice_host, lattice_device,
		domain_constants_host, domain_constants_device,
		domain_arrays_host, domain_arrays_device,
		output_controller_host, output_controller_device,
		times, project);*/
	int z_len = 1;
	#if DIM > 2
		z_len = domain_constants_host->length[2];
	#endif
	CGNSOutputHandler tmp("LBM-C Results.cgns",domain_constants_host->length[0],domain_constants_host->length[1],z_len);
	output_handler = tmp;
}



// COPIES f_i DATA FROM DEVICE TO HOST AND COMPUTERS MACROSCOPIC VALUES ON HOST, THIS DATA
// IS THEN WRITTEN TO THE OUTPUT FILE
//
// Note:	A computationally more efficient implementation would compute macroscopic
//			value's on the gpu and then just copy that data, this would however consume
//			more memory
void output_macros(int time)
{
	int domain_size = domain_constants_host->length[0]*domain_constants_host->length[1];
	#if DIM > 2
		domain_size = domain_size*domain_constants_host->length[2];
	#endif

	Lattice lattice_tmp;

	cudasafe(cudaMemcpy(&lattice_tmp, lattice_device, sizeof(Lattice),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
	
	double *u_tmp[DIM];
	cudasafe(cudaMemcpy(u_tmp, lattice_tmp.u, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	for(int d=0;d<DIM;d++)
	{
		cudasafe(cudaMemcpy(lattice_host->u[d], u_tmp[d], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
	}

	double *rho_tmp;
	cudasafe(cudaMemcpy(lattice_host->rho, lattice_tmp.rho, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	// Copy data from device to host
	//cudasafe(cudaMemcpy(lattice_host->rho, lattice_device->rho, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - rho");
	//cudasafe(cudaMemcpy(lattice_host->u[0], lattice_device->u[0], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - ux");
	//cudasafe(cudaMemcpy(lattice_host->u[1], lattice_device->u[1], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - uy");


	char *output_file;
	output_file = (char *)malloc(sizeof(char)*33);
	strcpy(output_file,"test3d.cgns");

	char **labels;
	double **data;

	int fields = 3;

	labels = (char **)malloc(fields * sizeof (char *));
	data = (double **)malloc(fields * sizeof(double));

	for(int i = 0; i<fields;i++)
	{
		labels[i] = (char *)malloc(STR_LENGTH*sizeof(char));
	}

	data[0] = lattice_host->rho;
	data[1] = lattice_host->u[0];
	data[2] = lattice_host->u[1];

	strcpy(labels[0],"Density");
	strcpy(labels[1],"VelocityX");
	strcpy(labels[2],"VelocityY");

	output_handler.append_solution_output(time,fields,data,labels);

	int i2d, i, j, k;
	i = 0;
	j = domain_constants_host->length[1]/2;
	#if DIM > 2
		k = domain_constants_host->length[2]/2;
	#else
		k = 0;
	#endif
	i2d = i+j*domain_constants_host->length[0]+k*domain_constants_host->length[0]*domain_constants_host->length[1];
	cout << endl << "time = " << time << "; rho = " << lattice_host->rho[i2d] << "; uX = " << lattice_host->u[0][i2d]<< "; uY = " << lattice_host->u[1][i2d] << "; resid = " << domain_constants_host->residual << endl;
}

// CONFIGURES THE KERNEL CONFIGURATION AND LAUNCHES KERNEL
void iterate(void)
{
	// GRID AND BLOCK DEFINITIONS CAN BE CALCULATED BEFORE ITERATE
	// DEFINE GRID AND BLOCK DIMS
	int3 threads;
	threads.x = (int)ceilf((float)domain_constants_host->length[0]/(float)NUM_THREADS_DIM_X);
	threads.y = (int)ceilf((float)domain_constants_host->length[1]/(float)NUM_THREADS_DIM_Y);
	threads.z = 1;

	int3 blocks;
	blocks.x = NUM_THREADS_DIM_X;
	blocks.y = NUM_THREADS_DIM_Y;
	blocks.z = 1;

	#if DIM >2
		threads.z = (int)ceilf((float)domain_constants_host->length[2]/(float)NUM_THREADS_DIM_Z);;
		blocks.z = NUM_THREADS_DIM_Z;
	#endif

	dim3 grid_dim = dim3(threads.x,threads.y,threads.z);
    dim3 block_dim = dim3(blocks.x,blocks.y,blocks.z);

	cudaThreadSynchronize();
	Check_CUDA_Error("Kernel \"iterate_bulk 1\" Execution Failed!");  
	// ITERATE ONCE
	iterate_kernel<<<grid_dim, block_dim>>>(lattice_device, domain_arrays_device, domain_constants_device, store_macros);
	cudaThreadSynchronize();
	Check_CUDA_Error("Kernel \"iterate_bulk 1\" Execution Failed!");  

	// SWAP CURR AND PREV LATTICE POINTERS READY FOR NEXT ITER
	swap_lattices();
}

void swap_lattices(void)
{
	cudasafe(cudaMemcpy(lattice_device_prototype, lattice_device, sizeof(Lattice),cudaMemcpyDeviceToHost),"Copy Data: Device Lattice Pointers From Device");
	double **tmp_1 = lattice_device_prototype->f_prev;
	double **tmp_2 = lattice_device_prototype->f_curr;
	lattice_device_prototype->f_curr = tmp_1;
	lattice_device_prototype->f_prev = tmp_2;
	cudasafe(cudaMemcpy(lattice_device, lattice_device_prototype, sizeof(Lattice),cudaMemcpyHostToDevice),"Copy Data: Device Lattice Pointers To Device");
}

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