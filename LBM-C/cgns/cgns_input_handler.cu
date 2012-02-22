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
	int index_file,index_base,index_zone,index_flow,index_field,index_coord;
	int icelldim, iphysdim, dim;

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
	CGNSInputHandler (char *);

	CGNSInputHandler ();

	template<class T>
	void read_field(T *data, char *label)
	{
		int num_arrays;

		open_file();

		cgns_error_check(cg_narrays(&num_arrays));

		for(int i = 0; i<num_arrays; i++)
		{
			// recover dataset labels from input file
			char *array_name;
			cg_array_info(i, array_name, NULL, NULL, NULL);

			//check for match and read or error
			if(*array_name == *label) 
			{
				cgns_error_check(cg_array_read(i,data));
				break;
			} else 
			{
				cout << endl << "Input Handler: " << label << " not found in file \"" << fname << "\"" << endl;
				exit(-1);
			}
		}

		close_file();

		cout << endl << "Input Handler: " << label << " loaded" << endl;
	}

};

CGNSInputHandler::CGNSInputHandler (char *input_filename) 
{
	fname = input_filename;

	#if DIM > 2
		icelldim=3;
		iphysdim=3;
		dim = 3;
	#else
		icelldim=2;
		iphysdim=2;
		dim = 2;
	#endif

	open_file();

}

CGNSInputHandler::CGNSInputHandler (){}

#endif