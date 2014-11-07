subroutine print_cmatrix(A,name)
complex(kind=8), dimension(:,:), intent(in) :: A
character(len=*), intent(in) :: name
character(Len=1) :: tab = char(9)

integer :: i,j, n, m
n = size(A,1)
m = size(A,2)
write(*,'(a,a)') name," = "
do i=1,n
    do j=1,m
       write(*,"("//FMT//",a,"//FMT//",a,a,$)") real(A(i,j),8), tab, aimag(A(i,j)), "i", tab
    enddo
    write(*,*) ""
enddo
write(*,*) ""
end subroutine print_cmatrix