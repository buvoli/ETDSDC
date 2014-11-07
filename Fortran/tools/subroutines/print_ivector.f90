subroutine print_ivector(A,name)
integer(kind=8),  intent(in) :: A(:)
character(len=*), intent(in) :: name
integer :: i, n
n = size(A)
write(*,'(a,a)') name," = "
do i=1,n
    write(*,"("//FMT//")") real(A(i),8)
enddo
write(*,*) ""
end subroutine print_ivector