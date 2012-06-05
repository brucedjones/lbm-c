#ifndef MACROS
#define MACROS

#define CGNS_TYPE CG_FILE_HDF5
//#define CGNS_TYPE CG_FILE_ADF

#define STR_LENGTH 64

#define D3Q15
//#define D2Q9

#ifdef D2Q9

	#define DIM 2
	#define Q 9

	#define LOAD_E(e) {																									    \
						e[0][0]= 0;e[0][1]= 1;e[0][2]= 0;e[0][3]=-1;e[0][4]= 0;e[0][5]= 1;e[0][6]=-1;e[0][7]=-1;e[0][8]= 1;	\
						e[1][0]= 0;e[1][1]= 0;e[1][2]= 1;e[1][3]= 0;e[1][4]=-1;e[1][5]= 1;e[1][6]= 1;e[1][7]=-1;e[1][8]=-1;	\
					  }
	#define LOAD_OMEGA(omega) {omega[0]=4.0/9.0;omega[1]=1.0/9.0;omega[2]=1.0/9.0;omega[3]=1.0/9.0;omega[4]=1.0/9.0;omega[5]=1.0/36.0;omega[6]=1.0/36.0;omega[7]=1.0/36.0;omega[8]=1.0/36.0;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=3;opp[2]=4;opp[3]=1;opp[4]=2;opp[5]=7;opp[6]=8;opp[7]=5;opp[8]=6;}

	#define NUM_THREADS_DIM_X 32
	#define NUM_THREADS_DIM_Y 16

#endif

#ifdef D3Q15

	#define DIM 3
	#define Q 15

	#define LOAD_E(e) {																																												\
						e[0][0]= 0;e[0][1]= 1;e[0][2]=-1;e[0][3]= 0;e[0][4]= 0;e[0][5]= 0;e[0][6]= 0;e[0][7]= 1;e[0][8]=-1;e[0][9]= 1;e[0][10]=-1;e[0][11]= 1;e[0][12]=-1;e[0][13]= 1;e[0][14]=-1;  \
						e[1][0]= 0;e[1][1]= 0;e[1][2]= 0;e[1][3]= 1;e[1][4]=-1;e[1][5]= 0;e[1][6]= 0;e[1][7]= 1;e[1][8]=-1;e[1][9]= 1;e[1][10]=-1;e[1][11]=-1;e[1][12]= 1;e[1][13]=-1;e[1][14]= 1;  \
						e[2][0]= 0;e[2][1]= 0;e[2][2]= 0;e[2][3]= 0;e[2][4]= 0;e[2][5]= 1;e[2][6]=-1;e[2][7]= 1;e[2][8]=-1;e[2][9]=-1;e[2][10]= 1;e[2][11]= 1;e[2][12]=-1;e[2][13]=-1;e[2][14]= 1;  \
					  }
	#define LOAD_OMEGA(omega) {omega[0]=2.0/9.0;omega[1]=1.0/9.0;omega[2]=1.0/9.0;omega[3]=1.0/9.0;omega[4]=1.0/9.0;omega[5]=1.0/9.0;omega[6]=1.0/9.0;omega[7]=1.0/72.0;omega[8]=1.0/72.0;omega[9]=1.0/72.0;omega[10]=1.0/72.0;omega[11]=1.0/72.0;omega[12]=1.0/72.0;omega[13]=1.0/72.0;omega[14]=1.0/72.0;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=2;opp[2]=1;opp[3]=4;opp[4]=3;opp[5]=6;opp[6]=5;opp[7]=8;opp[8]=7;opp[9]=10;opp[10]=9;opp[11]=12;opp[12]=11;opp[13]=14;opp[14]=13;}

	#define NUM_THREADS_DIM_X 32
	#define NUM_THREADS_DIM_Y 4
	#define NUM_THREADS_DIM_Z 4

#endif

#endif