===========================================================================
This directory contains modules which contain all information for defining
the linear and nonlinear operators for a given equation. Equation modules
also contains the numerical parameters for experiments. It is possible to
create your own module for solving different equations. In order for it to
function with existing code, a module must contain the following public 
variables and subroutines:


==== Public Variables/Subroutines Needed For SingleRun.f90 ====


1.	Np	INTEGER
		Number of spatial points.	

2.	tspan	DOUBLE array, dimension(2)
		Contains time integration bounds.

3.	y0	COMPLEX array, dimension(Np)
		Initial condition

4.	L	SUBROUTINE
		Returns PDE Linear Operator.
		
		Arguments (Diagonal Implementation)
		
			lambda	(output) COMPLEX*16 array, dimension(Np)
				Vector corresponding to diagonal entries of PDE 
				linear operator

		Arguments (Nondiagonal Implementation)

			lambda	(output) COMPLEX*16 array, dimension(Np,Np)
				Matrix corresponding to PDE linear operator

5.	N 	SUBROUTINE
		Evaluate PDE nonlinear operator N(t,\phi)

		Arguments
		
			t	(input)	DOUBLE
				time variable

			phi	(input) DOUBLE array, dimension(Np)
				solution phi at time t

			N_out	(output) Double array, dimension(NP)
				N_out = N(t,phi)

6.	init	SUBROUTINE
		Called before numerical experiments. Can be use to initializes 
		internal variables like initial condition. Init has no arguments.


==== Public Variables/Subroutines Needed For ExperimentRun.f90: ====


1.	All Variables/Subroutines for SingleRun.f90

2.	Fs	DOUBLE array, allocatable, dimension(:)
		Function Evaluations to include in method comparison. 
		I.E. Each method will be tested with F(i) function evaluations
		where i = 1...size(F)

3.	error_filter	SUBROUTINE 

		Arguments
		
			t	(input)	DOUBLE
				time variable

			phi	(input) DOUBLE array, dimension(Np)
				solution phi at time t

			N_out	(output) Double array, dimension(NP)
				N_out = N(t,phi)

4. 	reference_methods	LOCIGAL array, dimension(3)
				Specifies which method to include when computing reference solution 
				reference_methods(1) = true => use ETDSDC
				reference_methods(2) = true => use IMEXSDC
				reference_methods(3) = true => use ETDRK


5.	eqn_name	CHARACTER(len=*), parameter
			character array containing equation name. Results will be
			stored in the subdirectory result/eqn_name		