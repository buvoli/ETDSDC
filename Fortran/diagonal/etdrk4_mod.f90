! ============================ Module Description ===========================
! etdrk4_mod Implementation of ETDRK4 by Cox & Mathews. Contains subroutines:
!   1. etdrk4      - ETDRK4 timestepping code
!   2. initETDC    - initializes ETD coefficients
! ===========================================================================

module etdrk4_mod

    use tools_mod,  only: dp, tic, toc
    use phi_mod,    only: phi
    implicit none

    contains

    ! =======================================================================
    ! ETDRK4   Implements ETDRK4 Numerical Time Integrator
    !
    ! Arguments
    !
    !   L       (input) COMPLEX*16 array, dimension(N)
    !           vector cooresponding to PDE linear operator
    !
    !   N       (input) SUBROUTINE:
    !           cooresponds to PDE nonlinear operator. Must be of the form:
    !               subroutine N(t,y_in,N_out)
    !                   real(dp),    intent(in)  :: t
    !                   complex(dp), intent(in)  :: y_in(:)
    !                   complex(dp), intent(out) :: N_out(:)
    !               end subroutine N
    !
    !   y_in    (input) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   y_out   (output) COMPLEX*16 array, dimension (N)
    !           Initial condition.
    !
    !   times   (output) DOUBLE array, dimension(4)
    !           times(1:2) = [cputime,clocktime] for timestepping loop
    !           times(3:4) = [cputime,clocktime] for initializing coefficients,
    ! =======================================================================

    subroutine etdrk4(L,N,tspan,y_in,Nt,y_out,times)

        ! Arguments
        complex(dp), intent(in)     :: L(:)
        real(dp),    intent(in)     :: tspan(2)
        complex(dp), intent(in)     :: y_in(:)
        integer,         intent(in) :: Nt
        complex(dp), intent(out)    :: y_out(size(y_in))
        real(dp),    intent(out)    :: times(4)
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
        integer :: i,nL
        real(dp) :: h,t
        complex(dp), allocatable :: y(:)
        complex(dp), allocatable, dimension(:) :: E, E2, A123, b1, b23, b4, Nv, Na, Nb, Nc, a, b, c

        ! Initialize Variables
        nL = size(L)
        allocate(E(nL), E2(nL), A123(nL), b1(nL), b23(nL), b4(nL), Nv(nL), Na(nL), Nb(nL), Nc(nL), a(nL), b(nL), c(nL))
        y = y_in
        ! Set stepsize
        h = (tspan(2)-tspan(1))/(real(Nt,8))
        t = tspan(1)
        ! initialize ETD coefficients
        call tic()
        call initETDC(L,h,E, E2, A123, b1, b23, b4)
        call toc(times(3:4))
        ! timestepping loop
        call tic
        do i=1,Nt
            call N(t,y,Nv);
            a  = E2*y + 0.5_dp*A123*Nv;
            call N(t+h/2.0_dp, a, Na);
            b  = E2*y + 0.5_dp*A123*Na;
            call N(t+h/2.0_dp, b, Nb);
            c  = E2*a + A123*(Nb - 0.5*Nv);
            call N(t+h, c, Nc);
            y  =E*y + b1*Nv + b23*(Na + Nb) + b4*Nc;
            t = t + h;
        enddo
        call toc(times(1:2))
        y_out = y

    end subroutine etdrk4

    ! =======================================================================
    ! INITETDC   Initializes ETD coefficients for ETDRK4 scheme.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimension(N)
    !       vector/scaler cooresponding to PDE linear operator
    !
    !   h   (input) DOUBLE 
    !       timestep
    !
    !   E, E2, A, b1, b23, b4   (output) COMPLEX*16 array, dimension(N,N)
    !                           coefficients for ETDRK4
    ! =======================================================================

    subroutine initETDC(L,h,E, E2, A, b1, b23, b4)

        ! Arguments
        complex(dp), intent(in) :: L(:)
        real(dp),    intent(in) :: h
        complex(dp), dimension(:), intent(out) :: E, E2, A, b1, b23, b4

        ! Local Variables
        complex(dp), allocatable :: P(:,:)
        allocate(P(4,size(L)))

        ! Initialize Coefficients
        call phi(h/2.0_dp * L,1,P(1:2,:))
            E2  = p(1,:)
            A   = h*p(2,:)
        call phi(h * L,3,P)
            E   = P(1,:)
            b1  = h*(4.0_dp*P(4,:) - 3.0_dp*P(3,:) + P(2,:))
            b23 = h*(-4.0_dp*P(4,:) + 2.0_dp*P(3,:))
            b4  = h*(4.0_dp*P(4,:) - P(3,:))
        deallocate(P)

    end subroutine initETDC

end module etdrk4_mod