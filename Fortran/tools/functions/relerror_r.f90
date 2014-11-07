function relerror_r(A,B) result(E)
real(kind=8),  intent(in)  :: A(:)
real(kind=8),  intent(in)  :: B(:)
real(kind=8)  :: E
E = maxval(abs(A - B))/maxval(abs(B))
end function relerror_r