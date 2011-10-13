#ifndef MACROS
#define MACROS

#define Q 15
#define BLOCK_SIZE 64

#if Q == 15						
	#define LOAD_EX(ex) {ex[0]=0;ex[1]=1;ex[2]=-1;ex[3]=0;ex[4]=0;ex[5]=0;ex[6]=0;ex[7]=1;ex[8]=-1;ex[9]=1;ex[10]=-1;ex[11]=1;ex[12]=-1;ex[13]=1;ex[14]=-1;}
	#define LOAD_EY(ey) {ey[0]=0;ey[1]=0;ey[2]=0;ey[3]=1;ey[4]=-1;ey[5]=0;ey[6]=0;ey[7]=1;ey[8]=-1;ey[9]=1;ey[10]=-1;ey[11]=-1;ey[12]=1;ey[13]=-1;ey[14]=1;}
	#define LOAD_EZ(ez) {ez[0]=0;ez[1]=0;ez[2]=0;ez[3]=0;ez[4]=0;ez[5]=1;ez[6]=-1;ez[7]=1;ez[8]=-1;ez[9]=-1;ez[10]=1;ez[11]=1;ez[12]=-1;ez[13]=-1;ez[14]=1;}
	#define LOAD_OMEGA(omega) {omega[0]=2.f/9.f;omega[1]=1.f/9.f;omega[2]=1.f/9.f;omega[3]=1.f/9.f;omega[4]=1.f/9.f;omega[5]=1.f/9.f;omega[6]=1.f/9.f;omega[7]=1.f/72.f;omega[8]=1.f/72.f;omega[9]=1.f/72.f;omega[10]=1.f/72.f;omega[11]=1.f/72.f;omega[12]=1.f/72.f;omega[13]=1.f/72.f;omega[14]=1.f/72.f;}
	#define LOAD_OPP(opp) {opp[0]=0;opp[1]=2;opp[2]=1;opp[3]=4;opp[4]=3;opp[5]=6;opp[6]=5;opp[7]=8;opp[8]=7;opp[9]=10;opp[10]=9;opp[11]=12;opp[12]=11;opp[13]=14;opp[14]=13;}
#endif

#endif