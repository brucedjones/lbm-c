#ifndef CGNS_OUTPUT_HANDLER
#define CGNS_OUTPUT_HANDLER

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

#define STR_LENGTH 31

class ModelBuilder
{
	int length[DIM];
	// DEVICE VARIABLE DECLARATION
	Lattice *lattice_d;
	DomainArray *domain_arrays_d;
	DomainConstant *domain_constants_d;
	double **f_1_d, **f_2_d, *rho_d, **u_d, *boundary_value_d, *geometry_d, **force_d; 
	int *boundary_type_d;

	// HOST VARIABLE DECLARATION
	Timing time_t;
	Lattice *lattice_h, *lattice_d_prototype;
	DomainArray *domain_arrays_h;
	DomainConstant *domain_constants_h;
	double **f_h, *rho_h, **u_h, *boundary_value_h, *geometry_h, **force_h;
	int *boundary_type_h;

	// SCALAR DECLARATION (PLATFORM AGNOSTIC)
	double tau, residual;
	double tolerance;
	int domain_size, maxT, saveT, steadyT, collision_type;

	// CONFIG FLAGS
	bool zhou_he;
	bool forcing;
	bool is2D;
	
	void memory_allocator()
	{
		int domain_data_size;
		domain_data_size = domain_size*sizeof(double);

		// Allocate container structures
		combi_malloc<Lattice>(lattice_h, lattice_d, sizeof(Lattice));
		combi_malloc<DomainArray>(domain_arrays_h, domain_arrays_d, sizeof(DomainArray));
		combi_malloc<DomainConstant>(domain_constants_h, domain_constants_d, sizeof(DomainConstant));

		// Allocate required arrays
		combi_malloc<double>(f_h, f_1_d, Q*sizeof(f_h));
		cudasafe(cudaMalloc((void **)&df_2_d,Q*sizeof(f_h)), "Model Builder: Device memory allocation failed!");
		for(int i=0;i<Q;i++)
		{
			combi_malloc<double>(f_h[i], f_1_d[i], Q*sizeof(f_h));
			cudasafe(cudaMalloc((void **)&df_2_d[i],Q*sizeof(f_h)), "Model Builder: Device memory allocation failed!");
		}
		combi_malloc<double>(rho_h, rho_d, domain_data_size);
		combi_malloc<double>(u_h, u_d, domain_data_size*DIM);
		combi_malloc<double>(geometry_h, geometry_d, domain_data_size);
		
		
		// Allocate option arrays
		if(forcing) combi_malloc<double>(force_h, force_d, domain_data_size*DIM);
		if(zhou_he)
		{
			combi_malloc<int>(boundary_type_h, boundary_type_d, domain_data_size);
			combi_malloc<double>(boundary_value_h, boundary_value_d, domain_data_size);
		}
	}

	void memory_assembler()
	{
		lattice_h->f_prev = f_h;
		lattice_h->f_curr = f_h;
		lattice_h->u = u_h;
		lattice_h->rho = rho_h;

		Lattice *lattice_d_tmp = (Lattice *)malloc(sizeof(Lattice));
		lattice_d_tmp->f_prev = f_d;
		lattice_d_tmp->f_curr = f_d;
		lattice_d_tmp->u = u_d;
		lattice_d_tmp->rho = rho_d;
		cudasafe(cudaMemcpy(lattice_d, lattice_d_tmp, sizeof(Lattice),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

		domain_arrays_h->boundary_type = boundary_type_h;
		domain_arrays_h->boundary_value = boundary_value_h;
		domain_arrays_h->geometry = geometry_h;
		domain_arrays_h->force = force_h;

		DomainArray *domain_array_d_tmp = (DomainArray *)malloc(sizeof(DomainArray));
		domain_arrays_d_tmp->boundary_type = boundary_type_d;
		domain_arrays_d_tmp->boundary_value = boundary_value_d;
		domain_arrays_d_tmp->geometry = geometry_d;
		domain_arrays_d_tmp->force = force_d;
		cudasafe(cudaMemcpy(domain_array_d, domain_array_d_tmp, sizeof(DomainArray),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}

	void memory_loader()
	{
		CGNSInputHandler input_handler(fname_cgns, is2D);
		input_handler.read_field(domain_arrays_h->geometry, 'Porosity');
		cudasafe(cudaMemcpy(domain_arrays_d->geometry, domain_arrays_h->geometry, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		
		if(forcing)
		{
			char force_labels[3][33];
			force_labels[1] = "ForceX";
			force_labels[2] = "ForceY";
			force_labels[3] = "ForceZ";

			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain_arrays_h->force[d], &force_labels[d]);
				cudasafe(cudaMemcpy(domain_arrays_d->force[d], domain_arrays_h->force[d], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
			}
		}

		if(zhou_he)
		{
			input_handler.read_field(domain_arrays_h->boundary_type, "BCType");
			cudasafe(cudaMemcpy(domain_arrays_d->boundary_type, domain_arrays_h->boundary_type, sizeof(int)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

			input_handler.read_field(domain_arrays_h->boundary_value, "BCValue");
			cudasafe(cudaMemcpy(domain_arrays_d->boundary_value, domain_arrays_h->boundary_value, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		}

		if(static_ic)
		{
			load_static_IC();
		}
	}

void load_static_IC()
{
	int index_i;
	double omega[Q];
	LOAD_OMEGA(omega);
	for(int i=0;i<Q;i++)
	{
		for(int index=0;index<(domain_size);index++)
		{
			lattice_h->f_curr[i][index] = 1.0*omega[i];
			cudasafe(cudaMemcpy(lattice_d->f_curr[i], lattice_h->f_curr[i], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
			cudasafe(cudaMemcpy(lattice_d->f_prev[i], lattice_h->f_curr[i], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		}
	}
}

	template<class T>
	void combi_malloc(T *host_pointer, T *device_pointer, size_t size)
	{
		host_pointer = (T *)malloc(size);
		cudasafe(cudaMalloc((void **)&device_pointer,size), "Model Builder: Device memory allocation failed!");
	}

	void read_config()
	{
	}

public:
	ModelBuilder (char *, int, int, int);

	ModelBuilder ();

	get_model(Lattice *lattice_host, Lattice *lattice_device, DomainConstant *domain_constants_host, DomainConstant *domain_constants_device, DomainArray *domain_arrays_host, DomainArray *domain_arrays_host, Timing time)
	{
		lattice_host = lattice_h;
		lattice_device = lattide_d;
		domain_constants_host = domain_constants_h;
		domain_constants_device = domain_constants_d;
		domain_arrays_host = domain_arrays_h;
		domain_arrays_device = domain_arrays_d;
		time = time_t;

	}

};

CGNSOutputHandler::CGNSOutputHandler (char *input_filename) 
{
	fname_config = input_filename;
}

CGNSOutputHandler::CGNSOutputHandler (){}

#endif