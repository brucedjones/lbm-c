#ifndef MACROS
#define MACROS

#define DIM 2
#define Q 9

#if Q == 9						
	#define LOAD_EX(ex) {ex[0]=0;ex[1]=1;ex[2]=0;ex[3]=-1;ex[4]=0;ex[5]=1;ex[6]=-1;ex[7]=-1;ex[8]=1;}
	#define LOAD_EY(ey) {ey[0]=0;ey[1]=0;ey[2]=1;ey[3]=0;ey[4]=-1;ey[5]=1;ey[6]=1;ey[7]=-1;ey[8]=-1;}
	#define LOAD_OMEGA(omega) {omega[0]=4.0/9.0;omega[1]=1.0/9.0;omega[2]=1.0/9.0;omega[3]=1.0/9.0;omega[4]=1.0/9.0;omega[5]=1.0/36.0;omega[6]=1.0/36.0;omega[7]=1.0/36.0;omega[8]=1.0/36.0;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=3;opp[2]=4;opp[3]=1;opp[4]=2;opp[5]=7;opp[6]=8;opp[7]=5;opp[8]=6;}
#endif

#if DIM == 2
	#define NUM_THREADS_DIM_X 23
	#define NUM_THREADS_DIM_Y 22
#endif

#endif