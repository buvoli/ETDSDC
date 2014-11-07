function fourierWaveNum(n,Lx) result(ks)
! Returns fourier wavenumbers for n modes on domain of size Lx
! === Parameters ===
! n  - (real) number of modes
! Lx - (real) domain size
! === Output ===
! ks - (array) array of wavenumbers
integer,      intent(in)  :: n
real(kind=8), intent(in)  :: Lx
real(kind=8), allocatable :: ks(:)
integer      :: i
real(kind=8) :: s
allocate(ks(n))
s = (2.d0*PI)/Lx; ! Fourier wave number scaling factor
do i=0,n-1
    if (i <= n/2) then
        ks(i+1) = s*i
    else
        ks(i+1) = s*(-n + i)
    end if
enddo

end function fourierWaveNum