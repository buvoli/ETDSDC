! ============================ Module Description ===========================
! imexsdc_mod Implementation of IMEXSDC method. Contains a 2 subroutines:
!   1. imexsdc  - IMEXSDC timestepping code
!   2. initIM   - initializes Integration Matrix
! ===========================================================================

module imexsdc_mod

    ! Module Parameters
    use omp_lib
    use tools_mod, only: dp, tic, toc, weights
    implicit none
    type imexsdc_settings
        integer :: m
        real(dp), allocatable :: tau(:)
    end type imexsdc_settings

    contains

    ! =======================================================================
    ! IMEXSDC   Implements ETDSDC Numerical Time Integrator
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array, dimension(n)
    !           vector cooresponding to PDE linear operator
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(kind=8),    intent(in)  :: t
    !                   complex(kind=8), intent(in)  :: y_in(:)
    !                   complex(kind=8), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   tspan   (input) DOUBLE array, dimension(2)
    !           contains left and right integration bounds
    !
    !   y_in    (input) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   options (input) ETDSDC_SETTINGS
    !           derived data type that contains one field:
    !               options%tau : DOUBLE (quadrature points)
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    !           times(3:4) = [cputime,clocktime] for initializing coefficients,
    ! =======================================================================


    subroutine imexsdc(L,N,tspan,y_in,Nt,options,y_out,times)

        ! Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(inout)      :: y_in(:)
        type(imexsdc_settings), intent(in)  :: options
        integer,         intent(in)         :: Nt
        real(dp),    intent(out)        :: times(4)
        complex(dp), intent(out)        :: y_out(size(y_in))
        ! define interface subroutine function N
        interface
            subroutine N(t,yh_in, N_out)
            import :: dp
            real(dp),    intent(in)  :: t
            complex(dp), intent(in)  :: yh_in(:)
            complex(dp), intent(out) :: N_out(:)
            end subroutine N
        end interface

        ! Local Variables
        integer :: i,j,k,nn,m,nL
        real(dp)    :: h,t, eta(size(options%tau)-1), IM(size(options%tau)-1,size(options%tau))
        complex(dp), allocatable :: phi(:,:), phi_n(:,:), IR(:,:),  D(:,:), n_new(:)
        real(dp), allocatable :: tau(:)
        tau = options%tau
        m   = options%m
        nn  = size(tau)
        nL  = size(L)
        h   = (tspan(2)-tspan(1))/(real(Nt,8))
        t   = tspan(1)
        eta = tau(2:nn) - tau(1:nn-1)
        allocate(phi(nn,nL), phi_n(nn, nL),IR(nn-1, nL),D(nn-1,nL),n_new(nL))

        ! Initialize Integration Matrix & Constants
        call tic()
        call initIM(h,tau,IM)
        do i=1,nn-1
            do j=1,nL
                D(i,j) = (1.d0 - h*(eta(i))*L(j))
            enddo
        enddo
        call toc(times(3:4))

        ! timestepping loop
        call tic()
        phi(1,:) = y_in
        do i=1,Nt
            ! Euler step
            do j=1,nn-1
                call N(t + h*tau(j),phi(j,:),phi_n(j,:))
                phi(j+1,:) = (phi(j,:) + h*eta(j)*phi_n(j,:))/D(j,:)
            enddo
            call N(t + h*tau(nn),phi(nn,:),phi_n(nn,:))
            ! Correction Sweeps
            do k=1,m
                !$omp parallel do
                do j=1,nL
                    IR(:,j) = matmul(IM,L(j)*phi(:,j)+phi_n(:,j))
                enddo
                phi(2,:) = (phi(1,:) - h*eta(1)*L*phi(2,:) + IR(1,:))/D(1,:)
                do j=2,nn-1
                    call N(t+h*tau(j),phi(j,:),n_new)
                    phi(j+1,:) = (phi(j,:) - h*eta(j)*L*phi(j+1,:) + h*eta(j)*(n_new - phi_n(j,:)) + IR(j,:))/D(j,:)
                    phi_n(j,:) = n_new;
                enddo
                if(k /= m) then
                    call N(t + h*tau(nn),phi(nn,:),phi_n(nn,:))
                endif
            enddo
            phi(1,:) = phi(nn,:)
            t = t + h
        enddo
        call toc(times(1:2))
        y_out = phi(nn,:);
        
    end subroutine imexsdc

    ! =======================================================================
    ! INITIM  Initializes IMEXSDC integration matrix using weights function
    !         by Fornberg.
    !
    ! Arguments
    !
    !   h   (input) DOUBLE
    !       timestep
    !
    !   tau (input) DOUBLE array, dimensions(n,1)
    !       normalized quadrature points
    !
    !   IM  (output) DOUBLE array, dimension(n-1,n)
    !       integration matrix
    ! =======================================================================

    subroutine initIM(h,tau,IM)

        ! Arguments
        real(dp), intent(in)  :: h
        real(dp), intent(in)  :: tau(:)
        real(dp), intent(out) :: IM(size(tau)-1,size(tau))

        ! Local Variables
        real(dp) :: eta(size(tau)-1)
        real(dp) :: q(size(tau))
        real(dp) :: A(size(tau),size(tau))
        real(dp) :: p0(size(tau))
        integer  :: i,j,n

        ! Form Integation Matrix
        n     = size(tau)
        eta   = tau(2:n) - tau(1:n-1)

        ! Compute phi functions for Lambda=0 (taylor terms)
        p0(1) = 1
        do j=2,n
            p0(j) = p0(j-1)/j
        enddo

        do i=1,n-1
            q = (tau - tau(i))/eta(i)                   ! scaled quadrature points
            call weights(0.d0,q,n-1,A)                  ! finite difference matrix
            IM(i,:)   = h*eta(i)*matmul(a,p0)           ! store ith row of IM matrix
        enddo
    end subroutine initIM

end module imexsdc_mod