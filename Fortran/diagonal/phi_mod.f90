! ============================ Module Description ===========================
! PHI_MOD: Initilizes Phi functions for scaler arguments. Contains:
!   1. phi - initilizes phi functions using contour integral & recursion
!            relation.
! ===========================================================================

module phi_mod

    ! Module Parameters
    use tools_mod, only: dp, norm, PI, II
    implicit none

    ! Cauchy Integral Settings
    real(dp),   parameter :: c_tol = 1.0_dp         ! Smallest Lambda for contour integral
    real(dp),   parameter :: c_R   = 2.0_dp*c_tol   ! Contour Radius
    integer,    parameter :: c_M   = 32             ! Number of Points

    contains

    ! =======================================================================
    ! PHI   Evaluates \phi_i(L) for i=0,...,n and scaler/vector L using
    !       Recursion relation and Cauchy Integral Formula.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array, dimensions(n)
    !       array of Lambda values cooresponding to PDE linear component
    !
    !   n   (input) INTEGER
    !       highest phi function to initialize.
    !
    !   P   (output) COMPLEX*16, dimensions(n,size(L))
    !       array where on exit P(i,j) = \phi_{i-1}(L(j))
    ! =======================================================================

    subroutine phi(L,n,P)
        ! Arguments
        complex(dp), intent(in)  :: L(:)
        integer,     intent(in)  :: n
        complex(dp), intent(out) :: P(n+1,size(L))
        ! Local Variables
        integer     :: i,j,k,nL
        real(dp)    :: f
        complex(dp) :: z(c_M)
        complex(dp) :: Li,Lzi,pp

        nL = size(L)
        P = 0;
        ! Set contour points
        do i=1,c_M
            z(i) = c_R * exp(2.0_dp*PI*II*(i-1.0_dp)/real(c_M,dp))
        enddo
        ! Compute Phi
        do i=1,nL
            Li = L(i)
            if(abs(Li) >= c_tol) then
                ! Direct Formula
                P(1,i) = exp(Li)
                f = 1.0_dp;
                do j=2,n+1
                    P(j,i) = (P(j-1,i) - 1.0_dp/f)/Li
                    f = f*(j-1)
                enddo
            else
                ! Cauchy Integral Formula
                do k=1,c_M
                    Lzi = Li + z(k)
                    pp = exp(Lzi)
                    P(1,i) = P(1,i) + pp/c_M
                    f = 1.0_dp;
                    do j=2,n+1
                        pp = (pp - 1.0_dp/f)/Lzi
                        P(j,i) = P(j,i) + pp/c_M;
                        f = f*(j-1)
                    enddo
                enddo
                ! remove imaginary roundoff if L(i) is real
                if(aimag(Li) == 0.0_dp) then
                  P(:,i) = REALPART(P(:,i))
                endif
            end if
        end do
end subroutine phi

end module phi_mod