#ifndef INDEX_H
#define INDEX_H

void setup(char *data_file);
void cudasafe( cudaError_t error, char* message);
void Check_CUDA_Error(const char *message);
void output_macros(int time);
void iterate(void);
void swap_lattices(void);

double current_RMS(double *device_var, int var_size);
double error_RMS(double *device_var, int var_size);
void compute_residual(void);
void screen_mess(int iter, int coord[DIM]);

#endif