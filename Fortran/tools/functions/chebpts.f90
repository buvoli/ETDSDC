function chebpts(n,I_in) result(c)
! Returns chebyshev points on interval I
! === Parameters ===
! n - number of chebyshev points
! I - (array) 2x1 array containing inteval
! === Output ===
! c - (array) nx1 array of chebyshev points
integer, intent(in) :: n
real(kind=8), intent(in), optional :: I_in(2)
real(kind=8) :: c(n), m, a, I(2)
integer :: j
if(present(I_in)) then
    I = I_in
else
    I = (/ -1.d0, 1.d0 /)
endif

m = (I(2) + I(1))/2.d0
a = (I(2) - I(1))/2.d0
do j=1,n
    c(j) = -1.d0 * a * cos((j-1.d0)*PI/(n-1.d0)) + m
enddo
end function chebpts