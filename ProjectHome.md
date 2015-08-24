# Overview #

Current Version: [0.1](http://code.google.com/p/lbm-c/downloads/detail?name=LBM-C-0.1.zip)

LBM-C is a [lattice Boltzmann](http://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) 2D and 3D fluid flow solver implemented using nVidia's [CUDA](http://www.nvidia.com/object/cuda_home.html) platform. LBM-C is written in CUDA C and is licensed under GPL v2, you are invited to use, alter and improve upon LBM-C at your own discretion. The objectives of this open source code are:

  * To provide a straight-forward CUDA based implementation of the lattice Boltzmann method for modelling fluid flows
  * Encourage collaboration through the nature of open source code
  * Exercise the power of graphics hardware for the purposes of CFD (computational fluid dynamics)

If you wish to contribute to LBM-C a good place to start is the [LBM-C google group](http://groups.google.com/group/lbm-c).

# Features #
  * Models may be configured through use of a simple matlab interface
  * Input constants are specified through ASCII input file
  * Input data in CGNS format (compatible with most post-processors, though input data format is not cgns compliant)
  * Output data in CGNS format (compatible with most post-processors)

# Modelling Capabilities #
  * Arbitrary geometries may be considered either by use of full-way bounceback or a simplified case of [Nobel & Torczynski's](http://adsabs.harvard.edu/abs/1998IJMPC...9.1189N) immersed moving boundary method. The use of this boundary condition allows for simple coupling with the Discrete element method (not yet implemented), as was demonstrated by [Feng et. al](http://onlinelibrary.wiley.com/doi/10.1002/nme.2689/pdf)
  * Boundary pressures may be imposed upon edges using [Zhou & He's](http://pof.aip.org/resource/1/phfle6/v9/i6/p1591_s1?isAuthorized=no) pressure boundary condition
  * MRT or BGK collision dynamics
  * Smagorinsky turbulence (BGK only)
  * Body forces following [Guo's](http://pre.aps.org/abstract/PRE/v65/i4/e046308) work.
  * Operates on a D2Q9 or D3Q15 lattice

# Further Work #
  * Complete the implementation of Zhou & He's boundary condition to consider
    * Imposed pressure on convex corner nodes
    * Imposed velocity on edge's and corners
  * Multiphase coupling capability, where a suitable multiphase model would ideally include purely local computation
  * [Discrete element method](http://en.wikipedia.org/wiki/Discrete_element_method) coupling for particle laden flow

# Whats New #
**10/12/12**    New Release version 0.1, includes body forces, smagorinsky turbulence, and a matlab interface for constructing and post-processing models

