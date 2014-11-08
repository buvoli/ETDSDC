==================================================================================================

This Directory Contains the Fortran Code For The Numerical Experiments from
Buvoli T. "ETD Spectral Deferred Correction Methods", 2014

==================================================================================================

Directory Structure: 

	diagonal 	- Contains implementations of ETDSDC, IMEXSDC, and ETDRK4 for diagonal 
			  linear operator Lambda. diagonal/equations contains relevant equation
			  files.

	nondiagonal 	- Contains implementations of ETDSDC, IMEXSDC, and ETDRK4 for nondiagonal 
			  linear operator Lambda. nondiagonal/equations contains relevant equation 
			  files.

	tools 		- Contains functions and subroutines belonging to fortran module tools_mod. 
			  This module is used throughout the codebase by both the diagonal and
			  nondiagonal implementations.

	results		- Code outputs results to specified path located in this directory. Contains
			  two MATLAB scripts for generating figures.

==================================================================================================	

Compiling Instructions:	

	Makefiles are located in both the diagonal and nondiagonal directories. Each makefile has 
	the targets "experimentrun" and "singlerun".

		singlerun	- compiles program singlerun.f90 into executable singlerun which 
				  tests ETDSDC, IMEXSDC, and ETDRK4 schemes on one equation with
				  specified settings and prints results to terminal. 

		experimentrun	- compiles program experimentrun.f90 which reproduces one of the
				  four numerical experiments included in the paper. 

	Assuming Dependencies are satisfied (See makefile for details), typical usage might be:

		make singlerun
		./singlerun.exe

		make experimentrun
		./experimentrun.exe 	  

==================================================================================================

Solving Different Equations

	To change the equation being solved, modify line 9 of singlerun.f90 or experimentrun.f90
	and change equationname_mod.f90 to desired equation module. To generate new equation 
	modules see README.txt inside equation folder.
	
==================================================================================================

Figure Generation:

	Performance data for numerical schemes is collected in Fortran and saved in text files in
	the folder results. The MATLAB script loadData can be used to generate the figures in paper
	Note: this script requires the MATLAB implementation to generate solution figure.