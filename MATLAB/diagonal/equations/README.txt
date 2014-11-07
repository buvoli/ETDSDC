===========================================================================
This directory contains folders which contain all information for defining
the linear and nonlinear operators for a PDE as well as the numerical 
parameters for experiments. To add your own equation, each folder must 
contain the following files.

==== REQUIRED MATLAB FILES ====

1.	init.m	SCRIPT
		Must contain initialize the following parameters

		LF	vector
			Contains diagonal entries of PDE linear operator.
		
		NF	function
			Returns PDE nonlinear operator.

			Parameters
			t      - (scaler) time.
			phi    - (vector) PDE solution at time t.
			params - (struct) any additional parameters required.

		y0	vector
			Initial condition

		tspan	vector
			Contains time integration bounds.

==== OPTIONAL MATLAB FILES ====

1. 	L.m	FUNCTION
		L returns the PDE linear operator

		Parameters
		params - (struct) any additional parameters required.

2.	N.m	FUNCTION
		N evaluates the PDE nonlinear operator

		Parameters
		t      - (scaler) time.
		phi    - (vector) PDE solution at time t.
		params - (struct) any additional parameters required. 	

