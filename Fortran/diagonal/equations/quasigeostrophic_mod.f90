! ============================ Module Description ===========================
! quasigeostrophic_mod stores experiment parameters and subroutines for 
! evaluating Linear and Nonlinear Operators for the quasigeostrophic equation
! when solving in Fourier space (Diagonal Lambda).
!
! Contains 5 subroutines:
!   1. init         - initializes spatial grid, fft plans, experiment 
!                     parameters.
!   2. ic           - sets initial condition in physical space
!   3. L            - returns diagonal entries of linear operator Lambda
!   4. N            - subroutine for evaluating nonlinear operator
!   5. buildDiffMat - build relevant differentiation matrices
!   6. error_filter - defines error function used for numerical experiments
! ===========================================================================

module quasigeostrophic_mod

    ! Module Parameters
    use tools_mod, only: dp, PI, II, logspace, plinspace, fourierwavenum, meshgrid, relerror_c
    implicit none
    ! Numerical Parameters
    integer,    parameter :: Nx                   = 2**8                          ! Number of x spatial points (Must be Even)
    integer,    parameter :: Ny                   = 2**8                          ! Number of y spatial points (Must be Even)
    integer,    parameter :: Np                   = Nx * Ny                       ! Total Number of spatial points
    real(dp),   parameter :: tspan(2)             = [ 0.0_dp, 5.0_dp ]            ! Time integration window
    logical,    parameter :: reference_methods(3) = [ .true., .true., .false. ]   ! Methods for Reference Solution (ETDSDC,IMEXSDC,ETDRK)
    integer,    parameter :: num_tests            = 16                            ! Number of Numerical Tests
    real,       parameter :: smallest_F           = 1.5e3_dp                      ! Smallest Number of Function Evaluations
    real,       parameter :: largest_F            = 1.0e5_dp                      ! Maximum Number of Function Evaluations
    ! Storage Arrays
    real(dp), parameter :: Lx       = 2.0_dp * PI
    real(dp), parameter :: Ly       = 2.0_dp * PI
    real(dp), parameter :: epsilon  = 1.0e-2_dp
    real(dp), parameter :: v        = 1.0e-14_dp
    real(dp), parameter :: beta     = 10.0_dp
    character(len=*), parameter :: eqn_name = "quasigeostrophic"                   ! Used to save results to correct folder
    ! Storage Arrays
    complex(dp), dimension(:), allocatable   :: y0
    real(dp),    dimension(:,:), allocatable :: xs,ys
    complex(dp), dimension(:,:), allocatable :: y,yh,P_H, P_LAP, P_LAP_H, N_temp
    ! Diff Matrices
    complex(dp), dimension(:,:), allocatable :: DX, DY, LAP, ILAP
    ! FFT Parameters
    integer(dp) plan_backward,plan_forward
    include "fftw3.f90"
    ! Additional Settings
    real(dp), allocatable :: Fs(:) ! Function Counts to test (set in init function)

    contains

    subroutine init()
        integer     :: i
        real(dp)    :: k,s
        real(dp)    :: h = Lx/Nx
        real(dp), allocatable :: xsr(:), ysr(:)
        ! allocate arrays
        allocate(y(Nx,Ny), yh(Nx,Ny), xs(Nx,Ny), ys(Nx,Ny), P_H(Nx,Ny), P_LAP(Nx,Ny), P_LAP_H(Nx,Ny), N_temp(Nx,Ny))
        allocate(xsr(Nx),ysr(Ny),y0(Np))
        ! setup spatial domain
        xsr = plinspace(-Lx/2.0_dp,Lx/2.0_dp,Nx)
        ysr = plinspace(-Ly/2.0_dp,Ly/2.0_dp,Ny)
        call meshgrid(xsr,ysr,xs,ys)
        ! set FFT plans
        call dfftw_plan_dft_2d_ (plan_forward,  Nx, Ny, y, yh, FFTW_FORWARD,  FFTW_MEASURE)
        call dfftw_plan_dft_2d_ (plan_backward, Nx, Ny, yh, y, FFTW_BACKWARD, FFTW_MEASURE)
        ! Set Differentiation Matrices
        call buildDiffMat()
        ! Set Initial condition
        call ic()
        ! Set Function Evalutations for tests
        Fs = logspace(log(smallest_f)/log(10.0_dp),log(largest_F)/log(10.0_dp),num_tests)
    end subroutine

    ! Form Initial Condition
    subroutine ic()
        y = (1.0_dp/8.0_dp) * exp(-8.0_dp*(2.0_dp*ys**2.0_dp + 0.5_dp*xs**2.0_dp - PI/4.0_dp)**2.0_dp)
        call dfftw_execute(plan_forward)
        y0 = reshape(LAP * yh, (/ Nx*Ny /) )
    end subroutine ic

    subroutine L(lambda)
        complex(dp), dimension(Np), intent(out) :: lambda
        lambda = reshape(-1.0_dp*beta*DX*ILAP - epsilon - v*(DX**8 + DY**8), (/ Np /))
    end subroutine L

    ! Nonlinear Operator
    subroutine N(t,LAP_P_H_in, N_out)
        real(dp),    intent(in) :: t
        complex(dp), dimension(:), intent(in)  :: LAP_P_H_in !Fourier transform of Laplacian if Phi
        complex(dp), dimension(:), intent(out) :: N_out
        integer :: i,j

        yh(:,:) = reshape(LAP_P_H_in, [Nx, Ny]) !NOTE yh(:,:) = instead of yh = to avoid segfault with fftw?
        ! Fourier Transform of Phi
        P_H = ILAP * yh
        ! Laplacian of Phi
        call dfftw_execute(plan_backward)
        P_LAP = (1.0_dp/Np) * y
        ! Compute N1
        yh = -Dy * P_H
        call dfftw_execute(plan_backward)
        y  = (1.0_dp/Np) * y * P_LAP
        call dfftw_execute(plan_forward)
        N_temp  = -DX * yh
        ! Compute N2
        yh = Dx * P_H
        call dfftw_execute(plan_backward)
        y  = (1.0_dp/Np) * y * P_LAP
        call dfftw_execute(plan_forward)
        N_out  = reshape(N_temp - (DY * yh),[Np])
    end subroutine N

    ! Differention Matrices
    subroutine buildDiffMat()
        real(dp), allocatable :: k_x(:), k_y(:)
        integer :: i
        allocate(DX(NX,NY), DY(NX,NY), LAP(NX,NY), ILAP(NX,NY))
        ! Get Fourier wave numbers
        k_x = fourierWaveNum(Nx,Lx)
        k_y = fourierWaveNum(Ny,Ly)
        ! Form Discrete D/DX Diff matrix
        do i=1,nx
            DX(:,i) = II * k_x(i)
        enddo
        ! Form Discrete D/DY Diff matrix
        do i=1,nx
            DY(i,:) = II * k_y(i)
        enddo
        LAP  = DX**2 + DY**2        ! Laplacian operator
        ILAP = 1.0_dp/LAP           ! Inverse Laplacian
        ILAP(1,1) = (0.0_dp,0.0_dp) ! Set (0,0) mode to zero
    end subroutine buildDiffMat


    ! Error Filter (Inf norm in physical space)
    function error_filter(estimate,exact) result(error)

        complex(dp), intent(in)  :: estimate(Np)
        complex(dp), intent(in)  :: exact(Np)

        real(dp)  :: error
        complex(dp), allocatable :: a(:),b(:)

        allocate(a(Np),b(Np))

        ! IFFT of estimate
        yh(:,:) = (1.0_dp / Np) * reshape(estimate, [Nx, Ny])
        call dfftw_execute(plan_backward)
        a(:) = reshape(y,[Np])

        ! IFFT of exact
        yh(:,:) = (1.0_dp / Np) * reshape(exact, [Nx, Ny])
        call dfftw_execute(plan_backward)
        b(:) = reshape(y,[Np])

        error = relerror_c(a,b)
    end function error_filter

end module quasigeostrophic_mod