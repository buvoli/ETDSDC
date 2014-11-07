subroutine print_rvector(A,name)
real(kind=8),  intent(in) :: A(:)
character(len=*), intent(in) :: name
integer :: i, n
n = size(A)
write(*,'(a,a)') name," = "
do i=1,n
    write(*,"("//FMT//")") A(i)
enddo
write(*,*) ""
end subroutine print_rvector