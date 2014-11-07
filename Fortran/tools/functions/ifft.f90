function ifft(yh) result(y)
! FFT function
complex(kind=8), intent(in)  :: yh(:)
integer(kind=8) :: plan_backward
complex(kind=8), allocatable :: y(:)
allocate(y(size(yh)))
call dfftw_plan_dft_1d_ (plan_backward, size(yh), yh, y, FFTW_BACKWARD, FFTW_ESTIMATE)
call dfftw_execute(plan_backward)
y = y/real(size(y),8)
end function ifft