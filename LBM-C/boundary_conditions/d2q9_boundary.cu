#ifndef D2Q9_BOUNDARY
#define D2Q9_BOUNDARY

// Necessary includes
#include "../macros.cu"
#include "d2q9_boundary.cuh"

// These files are only included to remove squiggly red lines in VS2010
#include "../data_types.cuh"
#include "cuda_runtime.h"
#include "d2q9_zh_defs.cu"
#include "d2q9_sf_eq_defs.cuh"

__device__ __constant__ boundary_condition boundary_conditions[6] = { zh_pressure_x, zh_pressure_X, NULL, NULL, NULL, NULL,
															sf_eq_pressure_x,sf_eq_pressure_X,NULL, NULL, NULL, NULL};


#endif
