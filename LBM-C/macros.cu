#ifndef MACROS
#define MACROS

#define D2Q9

#ifdef D2Q9

	#define DIM 2
	#define Q 9

	#define LOAD_E(e) {e[0][0]=0;e[0][1]=1;e[0][2]=0;e[0][3]=-1;e[0][4]=0;e[0][5]=1;e[0][6]=-1;e[0][7]=-1;e[0][8]=1,e[1][0]=0;e[1][1]=0;e[1][2]=1;e[1][3]=0;e[1][4]=-1;e[1][5]=1;e[1][6]=1;e[1][7]=-1;e[1][8]=-1;}
	#define LOAD_OMEGA(omega) {omega[0]=4.0/9.0;omega[1]=1.0/9.0;omega[2]=1.0/9.0;omega[3]=1.0/9.0;omega[4]=1.0/9.0;omega[5]=1.0/36.0;omega[6]=1.0/36.0;omega[7]=1.0/36.0;omega[8]=1.0/36.0;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=3;opp[2]=4;opp[3]=1;opp[4]=2;opp[5]=7;opp[6]=8;opp[7]=5;opp[8]=6;}

	#define NUM_THREADS_DIM_X 20
	#define NUM_THREADS_DIM_Y 25

#endif

#ifdef D3Q15

	#define DIM 2
	#define Q 9

	#define LOAD_E(e) {e[0][0]=0;e[0][1]=1;e[0][2]=0;e[0][3]=-1;e[0][4]=0;e[0][5]=1;e[0][6]=-1;e[0][7]=-1;e[0][8]=1,
						e[1][0]=0;e[1][1]=0;e[1][2]=1;e[1][3]=0;e[1][4]=-1;e[1][5]=1;e[1][6]=1;e[1][7]=-1;e[1][8]=-1,
						e[2][0]=0;e[2][1]=0;e[2][2]=1;e[2][3]=0;e[2][4]=-1;e[2][5]=1;e[2][6]=1;e[2][7]=-1;e[2][8]=-1;}
	#define LOAD_OMEGA(omega) {omega[0]=4.0/9.0;omega[1]=1.0/9.0;omega[2]=1.0/9.0;omega[3]=1.0/9.0;omega[4]=1.0/9.0;omega[5]=1.0/36.0;omega[6]=1.0/36.0;omega[7]=1.0/36.0;omega[8]=1.0/36.0;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=3;opp[2]=4;opp[3]=1;opp[4]=2;opp[5]=7;opp[6]=8;opp[7]=5;opp[8]=6;}

	#define NUM_THREADS_DIM_X 20
	#define NUM_THREADS_DIM_Y 25

#endif

#endif