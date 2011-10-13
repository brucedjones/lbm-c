#ifndef INDEX_H
#define INDEX_H

void allocate_memory_host(void);
void allocate_memory_device(void);
void load_and_assemble_data(void);
void load_static_IC(void);
void setup(void);
void cudasafe( cudaError_t error, char* message);
void Check_CUDA_Error(const char *message);
void output_macros(int time);
void iterate(void);

#endif