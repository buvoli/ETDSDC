function fft2(y) result(yh)
! FFT2 function
complex(kind=8), intent(in)  :: y(:,:)
integer(kind=8) :: plan_forward
complex(kind=8), allocatable :: yh(:,:)
allocate(yh(size(y,1),size(y,2)))
call dfftw_plan_dft_2d_ (plan_forward,  size(y,1), size(y,2), y, yh, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute(plan_forward)
end function fft2