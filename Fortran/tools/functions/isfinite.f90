function isfinite(A) result(flag)
complex(kind=8), intent(in) :: A(:)
logical :: flag
real(kind=8), parameter :: INF = huge(1.d0)
real(kind=8) :: Ar,Ai
integer :: i,n

flag = .true.
n = size(A)
do i=1,n
    Ar = real(A(i))
    Ai = aimag(A(i))
    if((.not.(abs(Ar) < INF .and. abs(Ai) < INF)) .and. (isnan(Ar) .or. isnan(Ai))) then
        flag = .false.
        exit
    endif
enddo
end function isfinite