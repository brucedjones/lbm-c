#ifndef CGNS_OUTPUT_HANDLER
#define CGNS_OUTPUT_HANDLER


#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
/* cgnslib.h file must be located in directory specified by -I during compile: */


#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#define STR_LENGTH 31

class CGNSOutputHandler
{

	int ni_c,nj_c,nk_c;
	int ni_v, nj_v, nk_v;
	int vertex_num, cell_num;
	float *x, *y, *z;
	char *fname;

	vector<int> soltime;
	vector<string> solname;

	int *soltime_a;
	char *solname_a;

	// CGNS variables
	int index_file,index_base,index_zone,index_flow,index_field,index_coord;
	int icelldim, iphysdim;

	void create_file()
	{
		cg_set_file_type(CGNS_TYPE);
		cgns_error_check(cg_open(fname,CG_MODE_WRITE,&index_file));
		cgns_error_check(cg_close(index_file));
	}

	void open_file()
	{
/* open CGNS file for write */
	   cgns_error_check(cg_open(fname,CG_MODE_MODIFY,&index_file));
	}

	void close_file()
	{
/* close CGNS file */
		cgns_error_check(cg_close(index_file));
	}

	void write_base()
	{
		open_file();
		/* create base (user can give any name) */
		char basename[10];
	    strcpy(basename,"Base");
	    cgns_error_check(cg_base_write(index_file,basename,icelldim,iphysdim,&index_base));
	    close_file();

	    cout << endl << "Output Handler: CGNS output file " << fname << " succesfully initialised" << endl;
	}
	
	void write_grid()
	{
		int idx;

		cgsize_t isize2d[3][2], isize3d[3][3];

		allocate_grid_memory();
		open_file();
		
		for (int k=0; k < nk_v; k++)
		{
			for (int j=0; j < nj_v; j++)
			{
				for (int i=0; i < ni_v; i++)
				{
					if(nk_c==1)
					{
						idx = i+j*ni_v;
						x[idx]=i;
						y[idx]=j;
					} else {
						idx = i+j*ni_v+k*ni_v*nj_v;
						x[idx]=i;
						y[idx]=j;
						z[idx]=k;
					}
					
				}
			}
		}

		// define zone name (user can give any name) 
	    char zonename[33];
		strcpy(zonename,"LBM-C Output Data");
		
		if(nk_c==1) {
			// vertex size 
		    isize2d[0][0]=ni_v;
		    isize2d[0][1]=nj_v;
		    // cell size 
		    isize2d[1][0]=ni_c;
		    isize2d[1][1]=nj_c;
			// boundary vertex size (always zero for structured grids) 
		    isize2d[2][0]=0;
		    isize2d[2][1]=0;

			// create zone 
			cgns_error_check(cg_zone_write(index_file,index_base,zonename,*isize2d,Structured,&index_zone));
		} else {
			// vertex size 
			isize3d[0][0]=ni_v;
		    isize3d[0][1]=nj_v;
			isize3d[0][2]=nk_v;
			// cell size 
			isize3d[1][0]=ni_c;
		    isize3d[1][1]=nj_c;
		    isize3d[1][2]=nk_c;
			// boundary vertex size (always zero for structured grids) 
		    isize3d[2][0]=0;
		    isize3d[2][1]=0;
		    isize3d[2][2]=0;

			// create zone 
			cgns_error_check(cg_zone_write(index_file,index_base,zonename,*isize3d,Structured,&index_zone));
		}
		
	// write grid coordinates (user must use SIDS-standard names here) 
	   cgns_error_check(cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateX",
	       x,&index_coord));
	   cgns_error_check(cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateY",
	       y,&index_coord));
	   if(nk_c!=1){
	   cgns_error_check(cg_coord_write(index_file,index_base,index_zone,RealSingle,"CoordinateZ",
	       z,&index_coord));
	   }

		close_file();
		deallocate_grid_memory();

		cout << endl << "Output Handler: Grid succesfully written to " << fname << endl;
	}

	void allocate_grid_memory()
	{
		x = (float *) malloc(sizeof(float)*vertex_num);
		y = (float *) malloc(sizeof(float)*vertex_num);
		z = (float *) malloc(sizeof(float)*vertex_num);
	}

	void deallocate_grid_memory()
	{
		free(x);
		free(y);
		free(z);
	}

	void cgns_append_sol_field(double *field, char *name, int *index_field)
	{
		cgns_error_check(cg_field_write(index_file,index_base,index_zone,index_flow,RealDouble,name,field,index_field));
	}

	void write_iterative_data()
	{
	// TIME DEPENDENACE
	// create BaseIterativeData
		int nsteps=soltime.size();
		const cgsize_t nsteps_cg_tmp = soltime.size();
		const cgsize_t *nsteps_cg = &nsteps_cg_tmp;

		soltime_to_array();
		solname_to_array();

		cgns_error_check(cg_biter_write(index_file,index_base,"TimeIterValues",nsteps));
	// go to BaseIterativeData level and write time values
		cgns_error_check(cg_goto(index_file,index_base,"BaseIterativeData_t",1,"end"));
		cgns_error_check(cg_array_write("IterationValues",Integer,1,nsteps_cg,soltime_a));
	// create ZoneIterativeData
		cgns_error_check(cg_ziter_write(index_file,index_base,index_zone,"ZoneIterativeData_t"));
	// go to ZoneIterativeData level and give info telling which
	// flow solution corresponds with which time (solname(1) corresponds
	// with time(1), solname(2) with time(2), and solname(3) with time(3))
		//cgns_error_check(cg_goto(index_file,index_base,"ZoneIterativeData_t",1,"end"));
		cgns_error_check(cg_goto(index_file,index_base,"Zone_t",index_zone,"ZoneIterativeData_t",1,"end"));
		cgsize_t idata[2];
		idata[0] = STR_LENGTH+1;
		idata[1] = nsteps;
		cgns_error_check(cg_array_write("FlowSolutionPointers",Character,2,idata,solname_a));

	// add SimulationType
		cgns_error_check(cg_simulation_type_write(index_file,index_base,TimeAccurate));
		free_tmp_storage();
	}

	void soltime_to_array()
	{
		soltime_a = (int *)malloc(sizeof(int)*soltime.size());
		for(int i = 0; i<soltime.size(); i++)
		{
			soltime_a[i] = soltime.at(i);
		}
	}

	void solname_to_array()
	{
		char *name_tmp;
		string name_tmp_s;

		solname_a = (char *)malloc(soltime.size()*sizeof(char*)*(STR_LENGTH+1));
		for(int i = 0; i<soltime.size();i++)
		{
			//solname_a[i] = (char *)malloc(STR_LENGTH*sizeof(char));
			
			name_tmp_s = solname.at(i);

			name_tmp = (char *)malloc(sizeof(char)*name_tmp_s.size()+1);

			std::copy(name_tmp_s.begin(), name_tmp_s.end(), name_tmp);
			name_tmp[name_tmp_s.size()] = '\0';

			strcpy(&solname_a[i*(STR_LENGTH+1)],name_tmp);

		}
	}

	void free_tmp_storage()
	{
		/*for(int i = 0; i<soltime.size();i++)
		{
			free (solname_a[i]);
		}*/
		free (solname_a);
		free (soltime_a);
	}

	void cgns_error_check(int error_code)
	{
		if(error_code!=0)
		{
		const char *error_message = cg_get_error();
		cout << error_message << endl;
		getchar();
		cg_error_exit();
		}
	}

public:
	CGNSOutputHandler (char *, int, int, int);

	CGNSOutputHandler ();

	void append_solution_output(int iter, int num_fields, double **data, char **labels)
	{
		open_file();

		soltime.push_back(iter);

		stringstream node_name;//create a stringstream
		node_name << "Solution @ Time = " << iter;//add number to the stream
		string node_name_s = node_name.str();
		solname.push_back(node_name_s);

		char * node_name_c = new char[node_name_s.size() + 1];
		std::copy(node_name_s.begin(), node_name_s.end(), node_name_c);
		node_name_c[node_name_s.size()] = '\0';

		// create flow solution node
		cgns_error_check(cg_sol_write(index_file,index_base,index_zone,node_name_c, CellCenter,&index_flow));

		// write flow solution field
		for(int n = 0; n<num_fields;n++)
		{
			cgns_append_sol_field(data[n], labels[n], &index_field);
		}

		// write time dependant data
		write_iterative_data();

		close_file();
		cout << endl << "Output Handler: Solution @ Time = " << iter << " written" << endl;
	}

};

CGNSOutputHandler::CGNSOutputHandler (char *output_filename, int length_x, int length_y, int length_z) 
{
	fname = output_filename;// = output_filename;

	ni_c = length_x;
	nj_c = length_y;
	nk_c = length_z;

	cell_num = ni_c*nj_c*nk_c;

	ni_v = ni_c+1;
	nj_v = nj_c+1;
	nk_v = nk_c+1;

	

	if(nk_c==1)
	{
		icelldim=2;
		iphysdim=2;
		vertex_num = ni_v*nj_v;

	} else {
		icelldim=3;
		iphysdim=3;
		vertex_num = ni_v*nj_v*nk_v;
	}

	create_file();
	write_base();
	write_grid();
}

CGNSOutputHandler::CGNSOutputHandler (){}

#endif