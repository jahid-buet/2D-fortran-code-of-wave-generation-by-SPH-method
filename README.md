# 2D-fortran-code-of-wave-generation-by-SPH-method
This is a simple Fortran code of 2d wave generation by piston type wavemaker which is based on Smoothed Particles Hydrodynamics(SPH).
Source files of this fortran code is nearly same as my 1st repository "2D-Fortran-code-of-dam-break-using-SPH".Couple of modifications is done such as implementing moving piston type wavemaker,adding damping zone, determination of timestep size which was not implemented in my 1st repository etc.
To run the simulation one has to run the executable file "wave_generation_2d.exe". if you find any bug or problem in running the simulation u can notify me or u are free to change the source code situated in "Source files" folder. 
To modify the source code one has to recompile the source files in "Source" folder to create new executable file(.exe).The source files were compiled with ifort compiler.


# Definition of source files are summarized as follow:
- ComDense.f90: Compute rate of change of density of fluid particles using continuity equation and boundary particles using equation of state.

- ComGravity.f90: Compute acceleration due to gravity for all fluid particles in domain.

- Compress.f90: Compute pressure between fluid particles by using equation of state and boundary particles by using pressure extrapolated from fluid particles.

- ComPressGrad.f90: Compute acceleration due to pressure gradient for all particles(boundary+fluid) in domain.

- ComViscosity.f90: Compute acceleration due to vsicosity term for all fluid particles. Depending on  slip or no-slip boundary condition,viscosity for boundary particles are computed.
  Three different formulation of viscosity is implemented. Laminar+sps,artificial and orginal viscosity.User can choose any of one of the three viscosity formulation.

- Dam_break.f90: Main source file from where all subroutines are called  and user interface for input parametrs are defined .Also time integration is implemented in this source file.

- KernalGradientCorrection.f90: Compute inverse matrix related to formulation of kernal gradient correction which is used in only in viscosity(laminar+sps and artificial) computation.

- Time.f90:This module computes time step which is used in time integration  although in this  2d dam break simulation this routine is not used.

- Xsph.f90:This module computes Xsph to regularize the particles movement in SPH.it is not used in this simulation.

- geometry.f90:This module compute the required geometry(3 wall and 1 water column) according to the dimensions  given by user.

- initial.f90:This module is called from main source file before time initegration  and initialize variables(velocity,acceleration and pressure ) for all particles in domain.

- kernal.f90: This module computes kernel function needs for SPH simulation.Two types of kernal is implemented one is quintic spline kernal and other is wendland kernal.User can choose any one of the kernal for the simulation.

- kgf_matrix.f90:This module computes 3 by 3 inverse matrix for kernal gradient free SPH formulation  and only used in Orginal viscosity formulation.

- part.f90: This module define derived data types of particles and its various component.

- strain.f90:This module computes strain rate for fluid particles which is used in orginal viscosity formulation(direct discretiztion of 1st and 2nd order derivatives using SPH function and then kgf matrix are used to improve the accuracy of viscosity formulation).

- var.f90: This module is used to define initial density and compute constant mass of all particles in domain.

- vector.f90:This module contains different derived data types which is used by part.f90 module.
