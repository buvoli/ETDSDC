! Simplified calling sequence for lapack solver using existing LU decomposing (DGETRF)
function ZLUS(A,P,f) result(x)
    implicit none
    complex(kind=8), intent(in)    :: A(:,:)
    complex(kind=8), intent(in)    :: f(size(A,1))
    complex(kind=8), allocatable   :: x(:)
    integer, intent(in)         :: P(size(A,1))
    integer                     :: info,n

    n = size(A,1)
    allocate(x(n))
    x = f;
    call ZGETRS('N',n,1,A,n,P,x,n,INFO)  !http://www.netlib.no/netlib/lapack/complex16/zgetrs.f
end function ZLUS