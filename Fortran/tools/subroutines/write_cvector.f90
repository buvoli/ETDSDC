subroutine write_cvector(A,filename)

! Arguments
complex(kind=8), intent(in)  :: A(:)
character(len=*), intent(in) :: filename
character(Len=1)             :: tab = char(9)
! Local Variables
integer :: i, n, u
! write vector to file
n = size(A)
open(newunit=u,file=filename,status="replace")
do i=1,n
    write(u,"("//FMT//",a,"//FMT//",a)") real(A(i)), tab, aimag(A(i))
enddo
close(u)

end subroutine write_cvector