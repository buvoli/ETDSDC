! ============================ Module Description ===========================
! kdv_mod stores experiment parameters and subroutines for evaluating Linear
! and Nonlinear Operators for the KDV equation when solving in Physical space
! (Nondiagonal Lambda).
!
! Contains 6 subroutines:
!   1. init         - initializes spatial grid, fft plans, experiment 
!                     parameters.
!   2. ic           - sets initial condition in physical space
!   3. L            - returns diagonal entries of linear operator Lambda
!   4. N            - subroutine for evaluating nonlinear operator
!   5. DFTMatrix    - Constructs the forward or inverse DFT matrix
!   6. error_filter - defines error function used for numerical experiments
! ===========================================================================

module kdv_mod

    use tools_mod, only: dp, PI, II, logspace, plinspace, fourierwavenum, relerror_c
    implicit none
    ! Numerical Parameters
    integer,    parameter :: Np                   = 2**8                              ! Number of spatial points (Must be Even)
    real(dp),   parameter :: tspan(2)             = [ 0.d0, 3.6d0/PI ]                ! Time integration window
    logical,    parameter :: reference_methods(3) = [ .true., .true., .false. ]       ! Methods for Reference Solution (ETDSDC,IMEXSDC,ETDRK)
    integer,    parameter :: num_tests            = 16                                ! Number of Numerical Tests
    real,       parameter :: smallest_F           = 1.d3                              ! Smallest Number of Function Evaluations
    real,       parameter :: largest_F            = 2.d5                              ! Maximum Number of Function Evaluations
    ! Equation parameters
    real(dp), parameter     :: Lx       = 2.0_dp        ! Spatial Domain Size
    real(dp), parameter     :: delta    = 2.2e-2_dp
    character(len=*), parameter :: eqn_name = "kdv"     ! used to save results to correct folder
    ! Storage Arrays
    complex(dp), dimension(:), allocatable :: y,yh,ks,y0
    real(dp),    dimension(:), allocatable :: xs
    ! FFT Parameters
    integer(kind=8) plan_backward,plan_forward
    include "fftw3.f90"
    ! Additional Settings
    real(kind=8), allocatable :: Fs(:) ! Function Counts to test (set in init function)

    contains

    subroutine init()
        ! allocate arrays
        allocate(y(Np), yh(Np), ks(Np), xs(Np), Fs(num_tests))
        ! set x domain
        xs = plinspace(0.d0,Lx,Np)
        ! set FFT plans For L routine
        call dfftw_plan_dft_1d_ (plan_forward, Np, y, yh, FFTW_FORWARD, FFTW_MEASURE)
        call dfftw_plan_dft_1d_ (plan_backward, Np, yh, y, FFTW_BACKWARD, FFTW_MEASURE)
        ! set initial condition
        call ic()
        ! Set Fourier wave numbers
        ks = fourierWaveNum(Np,Lx)
        ks(Np/2+1) = 0;  ! zero out highest odd mode to insure real odd derivatives
        ! Set Function Evalutations for tests
        Fs = logspace(log(smallest_f)/log(10.d0),log(largest_F)/log(10.d0),num_tests)
    end subroutine

    ! Initial Condition
    subroutine ic()
        y0 = cos(PI * xs) ! NOTE: xs global var set in init
    end subroutine ic

    ! Linear Operator
    subroutine L(lambda)
        complex(kind=8), dimension(Np,Np), intent(out)  :: lambda
        complex(kind=8), allocatable, dimension(:,:)    :: DFT_F, DFT_B
        complex(kind=8), allocatable, dimension(:)      :: ls
        integer :: i,j
        ! Allocate DFT Matrices
        allocate(DFT_F(Np,Np), DFT_B(Np,Np), ls(Np))
        ! Set DFT Matrix
        call DFTMatrix(Np, .true.,  DFT_F)
        call DFTMatrix(Np, .false., DFT_B)
        ! Compute Fourier wave numbers
        ls = (-1.d0*(delta)**2)*((II*ks)**3)
        ! Form matrix product diag(ls) * DFT_F
        do i=1,Np
            DFT_F(:,i) = ls * DFT_F(:,i)
        enddo
        ! Form linear operator
        lambda = matmul(DFT_B,DFT_F)
        do i=1,Np
            do j=1,Np
                lambda(j,i) = complex(real(lambda(j,i),dp),0.0_dp)
            enddo
        enddo

    end subroutine L

    ! Nonlinear Operator
    subroutine N(t,y_in, N_out)
        real(kind=8),    intent(in) :: t
        complex(kind=8), dimension(:), intent(in)  :: y_in
        complex(kind=8), dimension(:), intent(out) :: N_out
        y = -0.5d0 * y_in**2
        call dfftw_execute(plan_forward)
        yh = II * ks * yh;
        call dfftw_execute(plan_backward)
        N_out = y/Np
    end subroutine N

    ! Constructs the forward or backward DFT matrix
    ! N       (integer) - number of points
    ! forward (boolean) - true=>forwards transform, false=>backwards transform
    ! A       (matrix)  - complex DFT Matrix
    subroutine DFTMatrix(N,forward,A)
        integer, intent(in) :: N
        logical, intent(in) :: forward
        complex(kind=8), intent(out) :: A(N,N)
        integer :: i,j
        complex(kind=8) :: omega

        if(forward) then
            omega = exp(-2.d0*Pi*II/real(N,8))
        else
            omega = exp(2.d0*Pi*II/real(N,8))
        endif
        ! Form Matrix
        do i=1,N
            do j=1,N
                A(i,j) = omega**((i-1)*(j-1))
            enddo
        enddo
        if(.not. forward) then
            A = 1.d0/real(N,8) * A
        endif
    end subroutine DFTMatrix

    ! Error Filter (Inf norm in physical space)
    function error_filter(estimate,exact) result(error)

        complex(dp), intent(in)  :: estimate(Np)
        complex(dp), intent(in)  :: exact(Np)

        real(dp)  :: error

        error = relerror_c(estimate,exact)
    end function error_filter

end module kdv_mod