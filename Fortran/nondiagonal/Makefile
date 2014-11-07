# ============================================================================================================
# Makefile for Fortran Implementation of Numerical Code in Buvoli T., "ETD Spectral Deferred Correction", 2014
# Complied using gfortran 4.9.0 and gnu make 3.81
#
# External Dependencies:
# 	BLAS 	(http://www.netlib.org/blas/)
# 	LAPACK 	(http://www.netlib.org/lapack/)
# 	FFTW 	(http://www.fftw.org)
# ============================================================================================================

# Include Appropriate Directories
VPATH = ../tools:equations

# Compiler
FC = gfortran

# Flags
FFLAGS  = -o3 -fexternal-blas
LDFLAGS = -L /usr/local/lib/ -lfftw3 -llapack -lblas

# Objects
EQNS = kdv_mod.o
METHODS = etdrk4_mod.o etdsdc_mod.o imexsdc_mod.o
PHI = phi_mod.o
TOOLS = tools_mod.o

PROGRAMS = singlerun experimentrun

all: $(PROGRAMS)

$(EQNS): $(TOOLS)
$(METHODS): $(TOOLS) $(PHI)

singlerun.o: $(TOOLS) $(METHODS) $(EQNS)
singlerun: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

experimentrun.o: $(TOOLS) $(METHODS) $(EQNS)
experimentrun: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

# ======================================================================
# General rules. Taken From
# (http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html)
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD