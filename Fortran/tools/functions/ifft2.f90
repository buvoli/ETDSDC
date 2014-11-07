function ifft2(yh) result(y)
! IFFT2 function
complex(kind=8), intent(in)     :: yh(:,:)
integer(kind=8)                 :: plan_backward
complex(kind=8), allocatable    :: y(:,:)
allocate(yh(size(y,1),size(y,2)))
call dfftw_plan_dft_2d_ (plan_backward, size(yh,1), size(yh,2), yh, y, FFTW_BACKWARD, FFTW_MEASURE)
call dfftw_execute(plan_backward)
y = y/real(size(y,1) * size(y,2),dp)
end function ifft2