#ifndef MODEL_BUILDER
#define MODEL_BUILDER

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
#include "infile_reader.cu"
#include "cgns/cgns_input_handler.cu"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda_util.cu"
using namespace std;

class ModelBuilder
{
	int length[DIM];
	// DEVICE VARIABLE DECLARATION
	OutputController *output_controller_d;
	Lattice *lattice_d;
	Domain *domain_d;
	DomainConstant *domain_constants_d;
	double **f_d, *rho_d, **u_d, *geometry_d, **force_d; 
	int *micro_bc_d;
	int *macro_bc_d;

	// HOST VARIABLE DECLARATION
	Timing *time_t;
	ProjectStrings *project_t;
	OutputController *output_controller_h;
	Lattice *lattice_h, *lattice_d_prototype;
	Domain *domain_h;
	DomainConstant *domain_constants_h;
	double **f_h, *rho_h, **u_h, *geometry_h, **force_h;
	int *micro_bc_h;
	int *macro_bc_h;

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
		//combi_malloc<Lattice>(&lattice_h, &lattice_d, sizeof(Lattice));
		//combi_malloc<Domain>(&domain_h, &domain_d, sizeof(Domain));
		//combi_malloc<DomainConstant>(&domain_constants_h, &domain_constants_d, sizeof(DomainConstant));
		//combi_malloc<OutputController>(&output_controller_h, &output_controller_d, sizeof(OutputController));
		//domain_constants_h = (DomainConstant *)malloc(sizeof(DomainConstant));
		//time_t = (Timing *)malloc(sizeof(Timing));
		//project_t = (ProjectStrings *)malloc(sizeof(ProjectStrings));
	}

	void constant_loader()
	{
		// LOAD CONSTANTS FROM FILE
		InfileReader infile_reader(fname_config, project_t, domain_constants_h, time_t, output_controller_h);
		
		// LOAD LATTICE CONSTANTS
		LOAD_E(domain_constants_h->e);
		LOAD_OMEGA(domain_constants_h->omega);
		LOAD_OPP(domain_constants_h->opp);
		LOAD_M(domain_constants_h->M);
		LOAD_M_INV(domain_constants_h->M_inv);
		for(int i =0;i<NUM_RESIDS;i++)
		{
			domain_constants_h->residual[i] = 1;
		}

		//transfer domain_constants to device (cant think of a better place to put this)
		//cudasafe(cudaMemcpy(domain_constants_d, domain_constants_h, sizeof(DomainConstant),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		cudasafe(cudaMemcpyToSymbol("domain_constants", domain_constants_h, sizeof(DomainConstant)),"Model Builder: Copy to device memory failed!");
		cudasafe(cudaMemcpy(output_controller_d, output_controller_h, sizeof(OutputController),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}

// Allocates memory for variables which have variable size due to problem geometry
	void variable_size_allocator()
	{
		domain_size = 1;
		for(int d = 0; d<DIM; d++)
		{
			domain_size = domain_size*domain_constants_h->length[d];
		}
		int domain_data_size;
		domain_data_size = domain_size*sizeof(double);

		// Allocate required arrays
		// PDFS
		double *f_tmp[Q];
		combi_malloc<double*>(&f_h, &f_d, sizeof(double*)*Q);
		for(int i=0;i<Q;i++)
		{
			combi_malloc<double>(&f_h[i], &f_tmp[i], domain_data_size);
		}
		cudasafe(cudaMemcpy(f_d,f_tmp,sizeof(double*)*Q,cudaMemcpyHostToDevice), "Model Builder: Device memory allocation failed!");

		// RHO
		combi_malloc<double>(&rho_h, &rho_d, domain_data_size);
		
		// VELOCITY
		double *u_tmp[DIM];
		combi_malloc<double*>(&u_h, &u_d, sizeof(double*)*DIM);
		for(int i=0;i<DIM;i++)
		{
			combi_malloc<double>(&u_h[i], &u_tmp[i], domain_data_size);
		}
		cudasafe(cudaMemcpy(u_d,u_tmp,sizeof(double*)*DIM, cudaMemcpyHostToDevice), "Model Builder: Device memory allocation failed!");

		// GEOMETRY
		combi_malloc<double>(&geometry_h, &geometry_d, domain_data_size);
		
		// ALLOCATE OPTION ARRAYS
		// FORCING
		if(domain_constants_h->forcing == true)
		{
			double *force_tmp[DIM];
			combi_malloc<double*>(&force_h, &force_d, sizeof(double*)*DIM);
			for(int i=0;i<DIM;i++)
			{
				combi_malloc<double>(&force_h[i], &force_tmp[i], domain_data_size);
			}
			cudasafe(cudaMemcpy(force_d,force_tmp,sizeof(double*)*DIM, cudaMemcpyHostToDevice), "Model Builder: Device memory allocation failed!");
		}

		// MICRO BC
		if(domain_constants_h->micro_bc == true)
		{
			combi_malloc<int>(&micro_bc_h, &micro_bc_d, domain_data_size);
		}

		// MACRO BC
		if(domain_constants_h->macro_bc == true)
		{
			combi_malloc<int>(&macro_bc_h, &macro_bc_d, domain_data_size);
		}
	}

	void variable_assembler()
	{
		lattice_h->f = f_h;


		Lattice *lattice_d_tmp = (Lattice *)malloc(sizeof(Lattice));
		lattice_d_tmp->f = f_d;
		cudasafe(cudaMemcpy(lattice_d, lattice_d_tmp, sizeof(Lattice),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

		domain_h->micro_bc = micro_bc_h;
		domain_h->macro_bc = macro_bc_h;
		domain_h->geometry = geometry_h;
		domain_h->force = force_h;
		domain_h->u = u_h;
		domain_h->rho = rho_h;

		Domain *domain_d_tmp = (Domain *)malloc(sizeof(Domain));
		domain_d_tmp->micro_bc = micro_bc_d;
		domain_d_tmp->macro_bc = macro_bc_d;
		domain_d_tmp->geometry = geometry_d;
		domain_d_tmp->force = force_d;
		domain_d_tmp->u = u_d;
		domain_d_tmp->rho = rho_d;
		cudasafe(cudaMemcpy(domain_d, domain_d_tmp, sizeof(Domain),cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}

	void variable_loader()
	{
		// LOAD GEOMETRY
		CGNSInputHandler input_handler(project_t->domain_fname, domain_constants_h->length);
		input_handler.read_field(domain_h->geometry, "Porosity");
		cudasafe(cudaMemcpy(geometry_d, geometry_h, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		
		// LOAD FORCES IF REQUIRED
		if(domain_constants_h->forcing == true)
		{
			char force_labels[3][33];
			strcpy(force_labels[0], "ForceX");
			strcpy(force_labels[1], "ForceY");
			strcpy(force_labels[2], "ForceZ");

			double *force_d_tmp[DIM];
			cudasafe(cudaMemcpy(force_d_tmp, force_d, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain_h->force[d], force_labels[d]);
				cudasafe(cudaMemcpy(force_d_tmp[d], force_h[d], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
			}
		}

		// LOAD MICRO BOUNDARY CONDITIONS IF REQUIRED
		if(domain_constants_h->micro_bc == true)
		{
			input_handler.read_field(domain_h->micro_bc, "MicroBC");
			cudasafe(cudaMemcpy(micro_bc_d, micro_bc_h, sizeof(int)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		}

		// LOAD MACRO BOUNDARY CONDITIONS IF REQUIRED
		if(domain_constants_h->macro_bc == true)
		{
			input_handler.read_field(domain_h->macro_bc, "MacroBC");
			cudasafe(cudaMemcpy(macro_bc_d, macro_bc_h, sizeof(int)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");

			char vel_labels[3][33];
			strcpy(vel_labels[0], "VelocityX");
			strcpy(vel_labels[1], "VelocityY");
			strcpy(vel_labels[2], "VelocityZ");

			double *u_d_tmp[DIM];
			cudasafe(cudaMemcpy(u_d_tmp, u_d, sizeof(double*)*DIM,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");
			for(int d=0;d<DIM;d++)
			{
				input_handler.read_field(domain_h->u[d], vel_labels[d]);
				cudasafe(cudaMemcpy(u_d_tmp[d], u_h[d], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
			}

			input_handler.read_field(domain_h->rho, "Rho");
			cudasafe(cudaMemcpy(rho_d, rho_h, sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
		}

		if(domain_constants_h->init_type == 0)
		{
			load_static_IC();
		}
	}

	void load_static_IC()
{
	double *f_d_tmp[Q];
	cudasafe(cudaMemcpy(f_d_tmp, f_d, sizeof(double*)*Q,cudaMemcpyDeviceToHost),"Model Builder: Copy from device memory failed!");

	double omega[Q];
	LOAD_OMEGA(omega);
	for(int i=0;i<Q;i++)
	{
		for(int index=0;index<(domain_size);index++)
		{
			lattice_h->f[i][index] = 1.0*omega[i];
		}
		cudasafe(cudaMemcpy(f_d_tmp[i], f_h[i], sizeof(double)*domain_size,cudaMemcpyHostToDevice),"Model Builder: Copy to device memory failed!");
	}
}

public:
	ModelBuilder (char *, Lattice*, Lattice*, DomainConstant*, DomainConstant*, Domain*, Domain*, OutputController*, OutputController*, Timing*, ProjectStrings*);

	ModelBuilder ();

};

ModelBuilder::ModelBuilder (char *input_filename, Lattice *lattice_host, Lattice *lattice_device, DomainConstant *domain_constants_host, DomainConstant *domain_constants_device, Domain *domain_host, Domain *domain_device, OutputController *output_controller_host, OutputController *output_controller_device, Timing *time, ProjectStrings *project) 
{
	lattice_h= lattice_host;
	lattice_d= lattice_device;
	domain_constants_h= domain_constants_host;
	domain_constants_d= domain_constants_device;
	domain_h= domain_host;
	domain_d= domain_device;
	output_controller_h= output_controller_host;
	output_controller_d = output_controller_device;
	time_t = time;
	project_t = project;

	fname_config = input_filename;
	constant_size_allocator();
	constant_loader();
	variable_size_allocator();
	variable_assembler();
	cout << "variable assembler complete" << endl;
	variable_loader();
	cout << "variable loader complete" << endl;
}

ModelBuilder::ModelBuilder (){}

#endif