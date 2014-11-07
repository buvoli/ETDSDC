module tools_mod

    implicit none
    ! FFT parameters
    include "misc/fftw3.f90"
    ! Double Precision
    integer, parameter :: dp=kind(0.d0)
    ! Print Output Format
    character(len=*), parameter  :: FMT = "ES23.15E3"
    ! Mathematical Constants
    real(dp),    parameter :: PI    = 4.0_dp*datan(1.0_dp)
    complex(dp), parameter :: II    = complex(0.0_dp,1.0_dp)
    ! Clock variables
    integer     :: clock_rate, start_clock, end_clock
    real(dp)    :: cpu_time_start, cpu_time_end

    contains

    ! FFT routines
    include "functions/fft.f90"
    include "functions/ifft.f90"

    ! Weights
    include "subroutines/weights.f90"

    ! File Output routines
    include "subroutines/write_cmatrix.f90"
    include "subroutines/write_rmatrix.f90"
    include "subroutines/write_cvector.f90"
    include "subroutines/write_rvector.f90"

    ! File Input routines
    include "subroutines/read_rvector.f90"

    ! Console Output routines
    include "subroutines/print_cmatrix.f90"
    include "subroutines/print_rmatrix.f90"
    include "subroutines/print_imatrix.f90"
    include "subroutines/print_cvector.f90"
    include "subroutines/print_rvector.f90"
    include "subroutines/print_ivector.f90"

    ! Matlab-like routines
    include "subroutines/tictoc.f90"
    include "subroutines/meshgrid.f90"

    include "functions/isfinite.f90"
    include "functions/linspace.f90"
    include "functions/logspace.f90"
    include "functions/norm.f90"

    ! Additional
    include "functions/plinspace.f90"
    include "functions/fourierwavnum.f90"
    include "functions/relerror_r.f90"
    include "functions/relerror_c.f90"
    include "functions/chebpts.f90"

    ! LAPACK & BLAS wrappers
    include "subroutines/ZLU.f90"
    include "subroutines/ZMM.f90"

    include "functions/ZLUS.f90"

end module tools_mod