#ifndef INDEX_H
#define INDEX_H

void setup(void);
void cudasafe( cudaError_t error, char* message);
void Check_CUDA_Error(const char *message);
void output_macros(int time);
void iterate(void);
void swap_lattices(void);

double current_RMS(double *device_var, int var_size);
double error_RMS(double *device_var, int var_size);

#endif