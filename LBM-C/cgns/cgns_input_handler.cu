#ifndef CGNS_INPUT_HANDLER
#define CGNS_INPUT_HANDLER

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

class CGNSInputHandler
{

	char *fname;

	// CGNS variables
	int index_file;
	int length[3];

	void open_file()
	{
		cgns_error_check(cg_open(fname,CG_MODE_READ,&index_file));
	}

	void close_file()
	{
		cgns_error_check(cg_close(index_file));
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
	CGNSInputHandler (char *, int [DIM]);

	CGNSInputHandler ();

	template<class T>
	void read_field(T *data, char *label)
	{
		int num_arrays;

		//unused pointers;
		DataType_t d_type;
		int d_dim;
		cgsize_t d_dim_vector;

		bool field_found = false;
		int i;
		char array_name[30];

		cgsize_t min[3], max[3];
		min[0] = 1;
		min[1] = 1;
		min[2] = 1;
		max[0] = length[0];
		max[1] = length[1];
		max[2] = length[2];

		open_file();

		cgns_error_check(cg_nfields(index_file, 1, 1, 1, &num_arrays));
		for(i = 1; i<num_arrays+1; i++)
		{
			cgns_error_check(cg_field_info(index_file, 1, 1, 1, i, &d_type, array_name));
			if(strcmp(array_name, label) == 0) 
			{
				field_found = true;
				cgns_error_check(cg_field_info(index_file, 1, 1, 1, i, &d_type, array_name));
				break;
			}
		}
		if(field_found==true)
		{
			cgns_error_check(cg_field_read(index_file, 1, 1, 1, label, d_type, min, max,data));
			cout << endl << "Input Handler: " << label << " loaded" << endl;
		} else {
			cout << endl << "Input Handler: " << label << " not found in file \"" << fname << "\"" << endl;
			exit(-1);
		}
		
		close_file();
		
	}

};

CGNSInputHandler::CGNSInputHandler (char *input_filename, int length_in[DIM]) 
{
	fname = input_filename;

	length[0] = length_in[0];
	length[1] = length_in[1];

	#if DIM > 2
		length[2] = length_in[2];
	#else
		length[2] = 1;
	#endif

	open_file();

}

CGNSInputHandler::CGNSInputHandler (){}

#endif