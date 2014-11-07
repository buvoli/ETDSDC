! Simplified calling sequence for lapack LU decomposition ZGETRF
! === Parameters ===
! A - NxN matrix to be factored (factorization will overwrite A)
! === Output ===
! P - permutation vector Nx1
subroutine ZLU(A,P)
    implicit none
    complex(kind=8), intent(inout)  :: A(:,:)
    integer,         intent(out)    :: P(size(A,1))
    integer :: info,n
    n = size(A,1)
    call ZGETRF(n,n,A,n,P,INFO) !http://www.netlib.no/netlib/lapack/complex16/zgetrf.f
end subroutine ZLU