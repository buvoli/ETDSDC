function fft(y) result(yh)
! FFT function
complex(kind=8), intent(in)  :: y(:)
integer(kind=8) :: plan_forward
complex(kind=8), allocatable :: yh(:)
allocate(yh(size(y)))
call dfftw_plan_dft_1d_ (plan_forward, size(y), y, yh, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan_forward)
end function fft