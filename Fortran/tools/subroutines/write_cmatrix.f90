subroutine write_cmatrix(A,filename)
complex(kind=8), intent(in)  :: A(:,:)
character(len=*), intent(in) :: filename
character(Len=1) :: tab = char(9)

integer :: i,j, n, m
n = size(A,1)
m = size(A,2)
open(unit=1,file=filename)

!write (1,"(I10,a,I10)") N,tab,M
do i=1,n
    do j=1,m
        write(1,"("//FMT//",a,"//FMT//",a,$)") real(A(i,j),8), tab, aimag(A(i,j)), tab
    enddo
    write(1,*) ""
enddo
write(1,*) ""

close(1)
end subroutine write_cmatrix