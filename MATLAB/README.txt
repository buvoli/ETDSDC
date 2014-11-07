==================================================================================================

	This Directory Contains the MATALB Code For The Numerical Methods from:
	Buvoli T. "ETD Spectral Deferred Correction Methods", 2014

==================================================================================================

Directory Structure: 

	diagonal 	- Contains implementations of ETDSDC, IMEXSDC, and ETDRK4 for diagonal 
			  linear operator Lambda. diagonal/equations contains relevant equation
			  files.

	nondiagonal 	- Contains implementations of ETDSDC, IMEXSDC, and ETDRK4 for nondiagonal 
			  linear operator Lambda. nondiagonal/equations contains relevant equation 
			  files.

	common 		- Contains function file for initializing phi functions, and a MATLAB 
			  implementation of Fornberg’s weight function.

	stability	- Contains MATLAB files for generating stability and accuracy  diagrams 
			  for ETDSDC and IMEXSDC

==================================================================================================	

Running Code:

	To solve a given differential equation simply run main.m in either the diagonal or 
	nondiagonal directories. To add your own equations, see README.txt inside equation folder.
	
==================================================================================================