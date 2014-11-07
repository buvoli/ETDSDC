subroutine print_cvector(A,name)
complex(kind=8),  intent(in) :: A(:)
character(len=*), intent(in) :: name
character(Len=1) :: tab = char(9)
integer :: i, n
n = size(A)
write(*,'(a,a)') name," = "
do i=1,n
    write(*,"("//FMT//",a,"//FMT//",a)") real(A(i),8), tab, aimag(A(i)), "i"
enddo
write(*,*) ""
end subroutine print_cvector