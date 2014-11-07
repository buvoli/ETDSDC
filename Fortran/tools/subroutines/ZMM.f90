! Simplified calling sequence for blas matrix multiplication ZGEMM on square matrices.
! Computes R = c * A*B where A,B matrices and c is a scaler.
! Arguments
!   A   (input) Complex*16, dimension(NxN) - Matrix A
!   B   (input) Complex*16, dimension(NxN) - Matrix B
!   c   (input) Complex*16 - scaler c
!   R   (output) Complex*16 - Matrix R = c*A*B

subroutine ZMM(A,B,c,R)
    implicit none
    complex(dp), intent(in)    :: A(:,:), B(:,:), c
    complex(dp), intent(out)   :: R(:,:)
    integer :: n
    n = size(A,1)
    ! For Details http://www.math.utah.edu/software/lapack/lapack-blas/zgemm.html
    call ZGEMM('N','N',n,n,n,c,A,n,B,n,0.0_dp,R,n)
end subroutine ZMM