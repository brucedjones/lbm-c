#ifndef CGNS_OUTPUT_HANDLER
#define CGNS_OUTPUT_HANDLER

#pragma comment(lib, "cgns/lib/cgns.lib")
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
/* cgnslib.h file must be located in directory specified by -I during compile: */
#include "cgns\include\cgnslib.h"
#include "cgns\cgns_output_handler.cuh"

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#define STR_LENGTH 31

class CGNSInputHandler
{

	char *fname;

	// CGNS variables
	int index_file,index_base,index_zone,index_flow,index_field,index_coord;
	int icelldim, iphysdim;

	void open_file()
	{
/* open CGNS file for read */
		index_base = 1;
		index_zone = 1;
		cgns_error_check(cg_open(fname,CG_MODE_READ,&index_file));
		cgns_error_check(cg_goto(index_file,index_base,"Zone_t",index_zone,"DiscreteData_t",1,"end"));
	}

	void close_file()
	{
/* close CGNS file */
		cgns_error_check(cg_close(index_file));
	}

	
	void cgns_append_sol_field(double *field, char *name, int *index_field)
	{
		cgns_error_check(cg_field_write(index_file,index_base,index_zone,index_flow,RealDouble,name,field,index_field));
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
	CGNSInputHandler (char *, bool);

	CGNSInputHandler ();

	void read_field(double *data, char *label)
	{
		open_file();
		//////////////////////neeed to do a cg_Array_read here
		cgns_error_check(cg_array_read(index_file,index_base,index_zone,label,RealDouble,irmin,irmax,data)
		cgns_error_check(cg_sol_write(index_file,index_base,index_zone,node_name_c, CellCenter,&index_flow));

		close_file();

		cout << endl << "Input Handler: " << label << " loaded" << endl;
	}

};

CGNSOutputHandler::CGNSInputHandler (char *input_filename, bool 2d) 
{
	fname = input_filename;

	if(2d)
	{
		icelldim=2;
		iphysdim=2;

	} else {
		icelldim=3;
		iphysdim=3;
	}

	open_file();

}

CGNSOutputHandler::CGNSOutputHandler (){}

#endif