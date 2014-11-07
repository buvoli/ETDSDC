function relerror_c(A,B) result(E)
complex(kind=8),  intent(in)  :: A(:)
complex(kind=8),  intent(in)  :: B(:)
real(kind=8)  :: E
E = maxval(abs(A - B))/maxval(abs(B))
end function relerror_c