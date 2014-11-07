! ============================ Module Description ===========================
! etdsdc_mod Implementation of ETDSDC method. Contains 2 subroutines:
!   1. etdsdc   - ETDSDC timestepping code
!   2. initW    - initializes ETD coefficients
! ===========================================================================

module etdsdc_mod

    ! Module Parameters
    use tools_mod,  only: dp, tic, toc, norm, weights
    use phi_mod,    only: phi

    implicit none
    ! Custom Settings Data Type
    type etdsdc_settings
        integer :: m
        real(dp), allocatable :: tau(:)
    end type etdsdc_settings

    contains

    ! =======================================================================
    ! ETDSDC    Implements ETDSDC Numerical Time Integrator
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array, dimension(M)
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

    subroutine etdsdc(L,N,tspan,y_in,Nt,options,y_out,times)

        !Arguments
        complex(dp), intent(in)         :: L(:)
        real(dp),    intent(in)         :: tspan(2)
        complex(dp), intent(in)         :: y_in(:)
        type(etdsdc_settings), intent(in)   :: options
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
        real(dp)    :: h,t
        complex(dp), allocatable :: W(:,:,:)
        complex(dp), allocatable :: P01(:,:,:)
        complex(dp), allocatable :: phi(:,:), phi_n(:,:), IW(:,:), n_new(:)
        real(dp),    allocatable :: tau(:)

        ! Initialize Variables
        tau = options%tau
        m   = options%m
        nn  = size(tau)
        nL  = size(L)
        h   = (tspan(2)-tspan(1))/(real(Nt,dp))
        t   = tspan(1)
        allocate(P01(2,nn-1,nL), W(nn-1,nn,nL), phi(nn,nL), phi_n(nn, nL), IW(nn, nL), n_new(nL))

        ! initialize ETD coefficients
        call tic()
        call initW(L,h,tau,W,P01)
        call toc(times(3:4))

        ! timestepping loop
        call tic()
        phi(1,:) = y_in
        do i=1,Nt
            ! Euler step
            do j=1,nn-1
                call N(t + h*tau(j),phi(j,:),phi_n(j,:))
                phi(j+1,:) = P01(1,j,:)*phi(j,:) + P01(2,j,:)*phi_n(j,:)
            enddo
            call N(t + h*tau(nn),phi(nn,:),phi_n(nn,:))
            ! Correction Sweeps
            do k=1,m
                do j=1,nL
                    IW(:,j) = matmul(W(:,:,j),phi_n(:,j))
                enddo
                phi(2,:) = P01(1,1,:)*phi(1,:) + IW(1,:)
                do j=2,nn-1
                    call N(t+h*tau(j),phi(j,:),n_new)
                    phi(j+1,:) = P01(1,j,:)*phi(j,:) + P01(2,j,:)*(n_new - phi_n(j,:)) + IW(j,:)
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
        y_out = phi(nn,:)

    end subroutine etdsdc

    ! =======================================================================
    ! INITW   Initializes ETDSDC W functions for matrix L using weights
    !         function by Fornberg.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimension(M)
    !       vector cooresponding to PDE linear operator
    !
    !   h   (input) DOUBLE 
    !       timestep
    !
    !   tau (input) DOUBLE array, dimensions(n,1)
    !       normalized quadrature points
    !
    !   W   (output) COMPLEX*16 array, dimension(size(tau)-1, size(tau), size(L))
    !       contains the ETD Coefficients so that W(:,:,j) is the ETD integration
    !       matrix cooresponding to lambda=L(j)
    !
    !   P01 (output) COMPLEX*16 array, dimension(2,size(tau)-1,size(L))
    !       contains phi_0 and phi_1 functions needed for ETD Euler method, where
    !       P01(:,:,1,k) = Exp(h_k L) and P01(:,:,2,k) = h*h_i*\phi_{1}(h_k L)
    ! =======================================================================

    subroutine initW(L,h,tau,W, P01)

        ! Arguments
        complex(dp), intent(in)  :: L(:)
        real(dp),    intent(in)  :: h
        real(dp),    intent(in)  :: tau(:)
        complex(dp), intent(out) :: W(size(tau)-1,size(tau),size(L))
        complex(dp), intent(out) :: P01(2,size(tau)-1,size(L))

        ! Local Variables
        real(dp) :: eta(size(tau)-1)
        real(dp) :: q(size(tau))
        real(dp) :: A(size(tau),size(tau))
        complex(dp), allocatable :: p(:,:)
        integer :: i,n

        ! Initialize variables
        n     = size(tau)
        eta   = tau(2:n) - tau(1:n-1)
        allocate(p(n+1,size(L)))

        ! Calculate W functions
        do i=1,n-1
            q = (tau - tau(i))/eta(i)                   ! scaled quadrature ppints
            call weights(0.d0,q,n-1,A)                  ! finite difference matrix
            call phi(eta(i)*h*L,n,p)                    ! phi functions 0 to n
            W(i,:,:)   = h*eta(i)*matmul(a,p(2:n+1,:))  ! store ith row of Integration matrix
            P01(1,i,:) = p(1,:)                         ! store exp(h_i L)
            P01(2,i,:) = h*eta(i)*p(2,:)                ! store h*eta(i)*phi_1(h_i L)
        enddo
        deallocate(p)

    end subroutine initW

end module etdsdc_mod