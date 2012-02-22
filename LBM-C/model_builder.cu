#ifndef MODEL_BUILDER
#define MODEL_BUILDER

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include "infile_reader.cu"
#include "cgns/cgns_input_handler.cu"
using namespace std;

#define STR_LENGTH 31

class ModelBuilder
{
	int length[DIM];
	// DEVICE VARIABLE DECLARATION
	OutputController *output_controller_d;
	Lattice *lattice_d;
	DomainArray *domain_arrays_d;
	DomainConstant *domain_constants_d;
	double **f_1_d, **f_2_d, *rho_d, **u_d, *boundary_value_d, *geometry_d, **force_d; 
	int *boundary_type_d;

	// HOST VARIABLE DECLARATION
	Timing *time_t;
	ProjectStrings *project_t;
	OutputController *output_controller_h;
	Lattice *lattice_h, *lattice_d_prototype;
	DomainArray *domain_arrays_h;
	DomainConstant *domain_constants_h;
	double **f_h, *rho_h, **u_h, *boundary_value_h, *geometry_h, **force_h;
	int *boundary_type_h;

	// SCALAR DECLARATION (PLATFORM AGNOSTIC)
	double tau, residual;
	double tolerance;
	int domain_size, maxT, saveT, steadyT, collision_type;

	// CONFIG FLAGS AND STRINGS
	char *fname_config;
	bool zhou_he;
	bool forcing;
	bool is2D;

// Allocates memory for variables which are constant in size
	void constant_size_allocator()
	{
		// Allocate container structures
		combi_malloc<Lattice>(lattice_h, lattice_d, sizeof(Lattice));
		combi_malloc<DomainArray>(domain_arrays_h, domain_arrays_d, sizeof(DomainArray));
		combi_malloc<DomainConstant>(domain_constants_h, domain_constants_d, sizeof(DomainConstant));
		combi_malloc<OutputController>(output_controller_h, output_controller_d, sizeof(OutputController));
		time_t = (Timing *)malloc(sizeof(Timing));
		project_t = (ProjectStrings *)malloc(sizeof(ProjectStrings));
	}

	void constant_loader()
	{
		InfileReader infile_reader(fname_config, project_t, domain_constants_h, time_t, output_controller_h);
		//transfer domain_constants to device (cant think of a better place to put this)
		cudasafe(cudaMemcpy(domain_constants_d, domain_constants_h, sizeof(DomainConstant),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		cudasafe(cudaMemcpy(output_controller_d, output_controller_h, sizeof(OutputController),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}

// Allocates memory for variables which have variable size due to problem geometry
	void variable_size_allocator()
	{
		int domain_data_size;
		domain_data_size = domain_size*sizeof(double);

		// Allocate required arrays
		combi_malloc<double*>(f_h, f_1_d, sizeof(f_h));
		cudasafe(cudaMalloc((void **)&f_2_d,Q*sizeof(f_h)), "Model Builder: Device memory allocation failed!");
		for(int i=0;i<Q;i++)
		{
			combi_malloc<double>(f_h[i], f_1_d[i], domain_data_size);
			cudasafe(cudaMalloc((void **)&f_2_d[i], domain_data_size), "Model Builder: Device memory allocation failed!");
		}
		combi_malloc<double>(rho_h, rho_d, domain_data_size);
		combi_malloc<double*>(u_h, u_d, sizeof(u_h));
		for(int i=0;i<DIM;i++)
		{
			combi_malloc<double>(u_h[i], u_d[i], domain_data_size);
		}
		combi_malloc<double>(geometry_h, geometry_d, domain_data_size);
		
		
		// Allocate option arrays
		if(domain_constants_h->forcing == true)
		{
			combi_malloc<double*>(force_h, force_d, sizeof(force_h));
			for(int i=0;i<DIM;i++)
			{
				combi_malloc<double>(force_h[i], force_d[i], domain_data_size);
			}
		}
		if(domain_constants_h->zhou_he == 1)
		{
			combi_malloc<int>(boundary_type_h, boundary_type_d, domain_data_size);
			combi_malloc<double>(boundary_value_h, boundary_value_d, domain_data_size);
		}
	}

	void variable_assembler()
	{
		lattice_h->f_prev = f_h;
		lattice_h->f_curr = f_h;
		lattice_h->u = u_h;
		lattice_h->rho = rho_h;

		Lattice *lattice_d_tmp = (Lattice *)malloc(sizeof(Lattice));
		lattice_d_tmp->f_prev = f_1_d;
		lattice_d_tmp->f_curr = f_2_d;
		lattice_d_tmp->u = u_d;
		lattice_d_tmp->rho = rho_d;
		cudasafe(cudaMemcpy(lattice_d, lattice_d_tmp, sizeof(Lattice),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

		domain_arrays_h->boundary_type = boundary_type_h;
		domain_arrays_h->boundary_value = boundary_value_h;
		domain_arrays_h->geometry = geometry_h;
		domain_arrays_h->force = force_h;

		DomainArray *domain_arrays_d_tmp = (DomainArray *)malloc(sizeof(DomainArray));
		domain_arrays_d_tmp->boundary_type = boundary_type_d;
		domain_arrays_d_tmp->boundary_value = boundary_value_d;
		domain_arrays_d_tmp->geometry = geometry_d;
		domain_arrays_d_tmp->force = force_d;
		cudasafe(cudaMemcpy(domain_arrays_d, domain_arrays_d_tmp, sizeof(DomainArray),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}

	void variable_loader()
	{
		// LOAD GEOMETRY
		CGNSInputHandler input_handler(project_t->domain_fname);
		input_handler.read_field(domain_arrays_h->geometry, "Porosity");
		cudasafe(cudaMemcpy(domain_arrays_d->geometry, domain_arrays_h->geometry, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		
		// LOAD FORCES IF REQUIRED
		if(domain_constants_h->forcing == true)
		{
			char force_labels[3][33];
			strcpy(force_labels[0], "ForceX");
			strcpy(force_labels[1], "ForceY");
			strcpy(force_labels[2], "ForceZ");

			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain_arrays_h->force[d], force_labels[d]);
				cudasafe(cudaMemcpy(domain_arrays_d->force[d], domain_arrays_h->force[d], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
			}
		}

		// LOAD ZHOU/HE VARIABLES IF REQUIRED
		if(domain_constants_h->zhou_he == 1)
		{
			input_handler.read_field(domain_arrays_h->boundary_type, "BCType");
			cudasafe(cudaMemcpy(domain_arrays_d->boundary_type, domain_arrays_h->boundary_type, sizeof(int)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

			input_handler.read_field(domain_arrays_h->boundary_value, "BCValue");
			cudasafe(cudaMemcpy(domain_arrays_d->boundary_value, domain_arrays_h->boundary_value, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		}

		if(domain_constants_h->init_type == 0)
		{
			load_static_IC();
		}
	}

	void load_static_IC()
{
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

public:
	ModelBuilder (char *);

	ModelBuilder ();

	void get_model(Lattice *lattice_host, Lattice *lattice_device, DomainConstant *domain_constants_host, DomainConstant *domain_constants_device, DomainArray *domain_arrays_host, DomainArray *domain_arrays_device, OutputController *output_controller_host, OutputController *output_controller_device, Timing *time, ProjectStrings *project)
	{
		lattice_host = lattice_h;
		lattice_device = lattice_d;
		domain_constants_host = domain_constants_h;
		domain_constants_device = domain_constants_d;
		domain_arrays_host = domain_arrays_h;
		domain_arrays_device = domain_arrays_d;
		output_controller_host = output_controller_h;
		output_controller_device = output_controller_d;
		time = time_t;
		project = project_t;
	}

};

ModelBuilder::ModelBuilder (char *input_filename) 
{
	fname_config = input_filename;
	constant_size_allocator();
	constant_loader();
	variable_size_allocator();
	variable_assembler();
	variable_loader();
}

ModelBuilder::ModelBuilder (){}

#endif