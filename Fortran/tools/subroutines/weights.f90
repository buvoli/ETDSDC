! =======================================================================
! WEIGHTS   Compute coefficients for finite difference approximation for
!           the derivatives 1 to m at point z assuming data is known at 
!           points in array x. Based on the program "weights" in 
!           B. Fornberg, "Calculation of weights in finite difference 
!           formulas", SIAM Review 40 (1998), pp. 685-691.
! Arguments
!
!   z   (input) DOUBLE
!       location where approximations are to be accurate
!
!   x   (input) DOUBLE Array
!       array containing interpolation points
!
!   m   (input) INTEGER
!       highest derivative for which weights are sought
!
!   W   (output) DOUBLE array, dimension(size(x),m+1)
!       matrix that gives weights at grid locations x for 
!       derivative of order j<=m are found in c(:,j)
! =======================================================================

subroutine weights(z,x,m,W)

    ! Arguments
    real(dp), intent(in)    :: z
    real(dp), intent(in)    :: x(:)
    integer,  intent(in)    :: m
    real(dp), intent(out)   :: W(size(x),m+1)

    ! Variable Declarations
    real(dp) :: c1, c2, c3, c4, c5
    integer  :: i,j,k,n,mn

    c1 = 1.0_dp
    c4 = x(1) - z
    W  = 0.0_dp
    W(1,1) = 1.0_dp

    n = size(x)
    do i=2,n
        mn = min(i,m+1)
        c2 = 1.0_dp
        c5 = c4
        c4 = x(i) - z
        do j=1,i-1
            c3 = x(i) - x(j)
            c2 = c2*c3;
            if(j == i-1) then
                do k=mn,2,-1
                    W(i,k) = c1*(real(k-1,dp)*W(i-1,k-1) - c5*W(i-1,k))/c2;
                enddo
                W(i,1) = -c1*c5*W(i-1,1)/c2;
            endif
            do k=mn,2,-1
                W(j,k) = (c4*W(j,k) - real(k-1,dp)*W(j,k-1))/c3;
            enddo
            W(j,1) = c4*W(j,1)/c3;
        enddo
        c1 = c2;
    enddo

end subroutine weights