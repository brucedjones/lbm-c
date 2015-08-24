# Instructions For Use #

  1. Install CUDA dev driver
  1. Install CUDA toolkit
  1. Run LBM-C.exe

CUDA files may be found at:

http://developer.nvidia.com/cuda-downloads

# Input File Format Specification #

```
X_Length Y_Length Relaxation_Time Output_Frequency Maximum_Timestep Initial_Condition
Boundary_Type Boundary_Value
Boundary_Type Boundary_Value
Boundary_Type Boundary_Value
...
```

Nodal boundary information is read in IJ order, that is I varies varies fastest, then J.

|Boundary\_Type|Description                                                  |
|:-------------|:------------------------------------------------------------|
|0             | Bounceback                                                  |
|1<=type<2     | Porous fluid node (1 = entirely fluid, 1.99 = Almost solid) |
|2             | Zhou/He prescribed pressure for x=0 and uy=0                |
|3             | Zhou/He prescribed pressure for x=lx-1 and uy=0             |