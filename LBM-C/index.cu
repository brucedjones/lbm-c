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
#include <thrust/transform_reduce.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/device_vector.h>

// DEVICE VARIABLE DECLARATION
Lattice *lattice_device;
Domain *domain_device;
DomainConstant *domain_constants_device;
OutputController *output_controller_device;

// HOST VARIABLE DECLARATION
Lattice *lattice_host, *lattice_device_prototype;
Domain *domain_host;
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

	// Get available memory on graphics card before allocation
	size_t freeMemory_before = 0;
	size_t totalMemory_before = 0;
	cudaMemGetInfo(&freeMemory_before, &totalMemory_before);
	
	// Initialise memory for LBM model
	setup(argv[1]);
	
	// Get available memory on graphics card after allocation
	size_t freeMemory_after = 0;
	size_t totalMemory_after = 0;
	cudaMemGetInfo(&freeMemory_after, &totalMemory_after);

	// Report program memory usage
	cout << "Total Device Memory:	 "<< totalMemory_after / 1024 / 1024 << "Mb" << endl;
	cout << "Total Availabe Memory:	 "<< freeMemory_before / 1024 / 1024 << "Mb" << endl;
	cout << "Memory Used:            "<< (freeMemory_before-freeMemory_after) / 1024 / 1024 << "Mb" << endl;

	// Report domain configuration
	printf("X-Length:		%d\n", domain_constants_host->length[0]);
	printf("Y-Length:		%d\n", domain_constants_host->length[1]);
	#if DIM > 2
		printf("Z-Length:		%d\n", domain_constants_host->length[2]);
	#endif
	printf("Relaxation Time (Tau):	%f\n", domain_constants_host->tau);
	printf("\nPress return to continue...");
	if (output_controller_host->interactive == true) getchar();

	// Get current clock cycle number
	clock_t t1=clock();

	int domain_size=1;
	int stop=0;
	for(int d = 0; d<DIM ;d++)
	{
		domain_size = domain_size*domain_constants_host->length[d];
	}

	for(int i = 1; i<times->max+1; i++)
	{
		if((times->plot>0 && i%times->plot == 0) ||
		   (times->steady_check>0 && i%times->steady_check) || 
		   (times->screen>0 && i%times->screen)) store_macros = true;

		iterate(i-1);

		if(times->plot>0 && i%times->plot == 0)
		{
			output_macros(i);
			store_macros = false;
		}

		if(times->screen>0 && i%times->screen == 0)
		{
			screen_mess(i,output_controller_host->screen_node);
			store_macros = false;
		}

		if(times->steady_check>0 && i%times->steady_check == 0)
		{
			compute_residual(i);
			
			for(int resid=0;resid<NUM_RESIDS;resid++)
			{
				if(domain_constants_host->residual[resid]<domain_constants_host->tolerance) stop += 1;
			}
			if(isIndeterminate(domain_constants_host->residual[i%NUM_RESIDS]))
			{
				output_macros(i);
				exit(1);
			} else if(stop==NUM_RESIDS)
			{
				output_macros(i);
				break;
			}
			stop = 0;
			store_macros = false;
		}
	}

	// Get current clock cycle number
	clock_t t2=clock();
	// Compare and report global execution time
	double cputime = ((double)t2-(double)t1)/(double)CLOCKS_PER_SEC;
	printf("\n\nTotal Run Time: %fs",cputime);
	printf("\nPress return to finish");
	if (output_controller_host->interactive == true) getchar();


}


// EXECUTES ALL ROUTINES REQUIRED FOR THE MODEL SET UP
void setup(char *data_file)
{
	// Set cuda device to use
	cudaSetDevice(0);
	cudaFuncSetCacheConfig(iterate_kernel, cudaFuncCachePreferL1);
	
	// Allocate container structures
	combi_malloc<Lattice>(&lattice_host, &lattice_device, sizeof(Lattice));
	combi_malloc<Domain>(&domain_host, &domain_device, sizeof(Domain));
	combi_malloc<DomainConstant>(&domain_constants_host, &domain_constants_device, sizeof(DomainConstant));
	combi_malloc<OutputController>(&output_controller_host, &output_controller_device, sizeof(OutputController));
	domain_constants_host = (DomainConstant *)malloc(sizeof(DomainConstant));
	times = (Timing *)malloc(sizeof(Timing));
	project = (ProjectStrings *)malloc(sizeof(ProjectStrings));
	lattice_device_prototype = (Lattice *)malloc(sizeof(Lattice));

	ModelBuilder tmpmb(data_file, lattice_host, lattice_device,
		domain_constants_host, domain_constants_device,
		domain_host, domain_device,
		output_controller_host, output_controller_device,
		times, project);
	model_builder = tmpmb;

	int z_len = 1;
	#if DIM > 2
		z_len = domain_constants_host->length[2];
	#endif
	CGNSOutputHandler tmp(project->output_fname,domain_constants_host->length[0],domain_constants_host->length[1],z_len);
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

	Domain domain_tmp;

	cudasafe(cudaMemcpy(&domain_tmp, domain_device, sizeof(Domain),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
	
	double *u_tmp[DIM];
	cudasafe(cudaMemcpy(u_tmp, domain_tmp.u, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	for(int d=0;d<DIM;d++)
	{
		cudasafe(cudaMemcpy(domain_host->u[d], u_tmp[d], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
	}

	double *rho_tmp;
	cudasafe(cudaMemcpy(domain_host->rho, domain_tmp.rho, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	// Copy data from device to host
	//cudasafe(cudaMemcpy(lattice_host->rho, lattice_device->rho, sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - rho");
	//cudasafe(cudaMemcpy(lattice_host->u[0], lattice_device->u[0], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - ux");
	//cudasafe(cudaMemcpy(lattice_host->u[1], lattice_device->u[1], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - uy");

	int num_fields = 0;
	if (output_controller_host->u[0] == true) num_fields++;
	if (output_controller_host->u[1] == true) num_fields++;
#if DIM > 2
	if (output_controller_host->u[2] == true) num_fields++;
#endif
	if (output_controller_host->rho == true) num_fields++;

	char **labels;
	double **data;

	labels = (char **)malloc(num_fields * sizeof (char *));
	data = (double **)malloc(num_fields * sizeof(double));

	for(int i = 0; i<num_fields;i++)
	{
		labels[i] = (char *)malloc(STR_LENGTH*sizeof(char));
	}

	int counter = 0;

	if (output_controller_host->u[0] == true)
	{
		data[counter] = domain_host->u[0];
		strcpy(labels[counter],"VelocityX");
		counter++;
	}

	if (output_controller_host->u[1] == true)
	{
		data[counter] = domain_host->u[1];
		strcpy(labels[counter],"VelocityY");
		counter++;
	}
#if DIM > 2
	if (output_controller_host->u[2] == true)
	{
		data[counter] = domain_host->u[2];
		strcpy(labels[counter],"VelocityZ");
		counter++;
	}
#endif	
	if (output_controller_host->rho == true)
	{
		data[counter] = domain_host->rho;
		strcpy(labels[counter],"Density");
		counter++;
	}

/*	data[0] = lattice_host->rho;
	data[1] = lattice_host->u[0];
	data[2] = lattice_host->u[1];

	strcpy(labels[0],"Density");
	strcpy(labels[1],"VelocityX");
	strcpy(labels[2],"VelocityY");*/

	output_handler.append_solution_output(time,num_fields,data,labels);
}

// CONFIGURES THE KERNEL CONFIGURATION AND LAUNCHES KERNEL
void iterate(int t)
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
	iterate_kernel<<<grid_dim, block_dim>>>(lattice_device, domain_device, store_macros,t);
	cudaThreadSynchronize();
	Check_CUDA_Error("Kernel \"iterate_bulk 1\" Execution Failed!");  
	// SWAP CURR AND PREV LATTICE POINTERS READY FOR NEXT ITER
	//swap_lattices();
}

#if DIM > 2
	struct energy
	{
	    template <typename Tuple>
	    __host__ __device__
	    void operator()(Tuple t)
	    {
	        thrust::get<4>(t) = 0.5*thrust::get<3>(t)*((thrust::get<0>(t)*thrust::get<0>(t)) + (thrust::get<1>(t)*thrust::get<1>(t)) + (thrust::get<2>(t)*thrust::get<2>(t)));
	    }
	};
#else
	struct energy
	{
	    template <typename Tuple>
	    __host__ __device__
	    void operator()(Tuple t)
	    {
	        thrust::get<3>(t) = 0.5*thrust::get<2>(t)*((thrust::get<0>(t)*thrust::get<0>(t)) + (thrust::get<1>(t)*thrust::get<1>(t)));
	}
	};
#endif


double current_RMS(double *device_var_u[DIM], double *device_var_rho, int var_size)
{
	double *result;
	cudasafe(cudaMalloc((void **)&result,sizeof(double)*var_size), "Model Builder: Device memory allocation failed!");

	// wrap raw pointer with a device_ptr for thrust compatibility
	thrust::device_ptr<double> dev_ptr_x(device_var_u[0]);
	thrust::device_ptr<double> dev_ptr_y(device_var_u[1]);
	#if DIM > 2
		thrust::device_ptr<double> dev_ptr_z(device_var_u[2]);
	#endif
	thrust::device_ptr<double> dev_ptr_rho(device_var_rho);
	thrust::device_ptr<double> dev_ptr_res(result);

	// apply the transformation
	#if DIM > 2
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(dev_ptr_x, dev_ptr_y, dev_ptr_z, dev_ptr_rho, dev_ptr_res)),
                     thrust::make_zip_iterator(thrust::make_tuple(dev_ptr_x+var_size, dev_ptr_y+var_size, dev_ptr_z+var_size, dev_ptr_rho+var_size, dev_ptr_res+var_size)),
                     energy());
	#else
		thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(dev_ptr_x, dev_ptr_y, dev_ptr_rho, dev_ptr_res)),
                     thrust::make_zip_iterator(thrust::make_tuple(dev_ptr_x+var_size, dev_ptr_y+var_size, dev_ptr_rho+var_size, dev_ptr_res+var_size)),
                     energy());
	#endif
	Check_CUDA_Error("Steady State Calculation Kernel Execution Failed!");  
    
	// Compute RMS value
	//double sum = thrust::reduce(dev_ptr_res, dev_ptr_res+var_size, (double) 0, thrust::plus<double>());
	//double curr_RMS = sqrt(sum/var_size);

	double curr_RMS = thrust::reduce(dev_ptr_res, dev_ptr_res+var_size, (double) 0, thrust::plus<double>());

	cudasafe(cudaFree(result),"Freeing Device Memory");

	return curr_RMS;
}

double prev_RMS = 0;

double error_RMS(double *device_var_u[DIM], double *device_var_rho, int var_size)
{
	double curr_RMS = current_RMS(device_var_u, device_var_rho, var_size);
	double tmp = ((abs(curr_RMS-prev_RMS)/times->steady_check))/curr_RMS;

	prev_RMS = curr_RMS;

	return tmp;
}

void compute_residual(int time)
{
	int domain_size = domain_constants_host->length[0]*domain_constants_host->length[1];
	#if DIM > 2
		domain_size = domain_size*domain_constants_host->length[2];
	#endif

	Domain domain_tmp;

	cudasafe(cudaMemcpy(&domain_tmp, domain_device, sizeof(Domain),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	double *u_tmp[DIM];
	cudasafe(cudaMemcpy(u_tmp, domain_tmp.u, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	//double *rho_tmp;
	//cudasafe(cudaMemcpy(rho_tmp, domain_tmp.rho, sizeof(double*),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	/*cudasafe(cudaMemcpy(domain_host->u[0], u_tmp[0], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");
	cudasafe(cudaMemcpy(domain_host->u[1], u_tmp[1], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");
	cudasafe(cudaMemcpy(domain_host->u[2], u_tmp[2], sizeof(double)*domain_size,cudaMemcpyDeviceToHost),"Copy Data: Output Data - u");*/

//	domain_constants_host->residual = error_RMS(u_tmp[0],u_tmp[1],u_tmp[2], rho_tmp,domain_size);
	domain_constants_host->residual[time%NUM_RESIDS] = error_RMS(u_tmp, domain_tmp.rho,domain_size);
}

void screen_mess(int iter, int coord[DIM])
{
	int idx = coord[0]+coord[1]*domain_constants_host->length[0];
	#if DIM > 2
		idx += coord[2]*domain_constants_host->length[0]*domain_constants_host->length[1];
	#endif

	double u[DIM],rho;
	Domain domain_tmp;

	cudasafe(cudaMemcpy(&domain_tmp, domain_device, sizeof(Domain),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
	
	double *u_tmp[DIM];
	cudasafe(cudaMemcpy(u_tmp, domain_tmp.u, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	for(int d=0;d<DIM;d++)
	{
		cudasafe(cudaMemcpy(&u[d], &u_tmp[d][idx], sizeof(double),cudaMemcpyDeviceToHost),"Model Builder: BLAHBLAHCopy from device memory failed!");
	}

	cudasafe(cudaMemcpy(&rho, &domain_tmp.rho[idx], sizeof(double),cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	cout << "time = " << iter << "; rho = " << rho << "; uX = " << u[0]<< "; uY = " << u[1] << "; ";
	#if DIM>2
		cout << "uZ = " << u[2] << "; ";
	#endif
	cout << "resid = " << domain_constants_host->residual[iter%NUM_RESIDS] << endl;
}

bool isIndeterminate(const double pV)
{
    return (pV != pV);
} 

#endif