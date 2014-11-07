subroutine write_rvector(A,filename)

! Arguments
real(dp), intent(in)          :: A(:)
character(len=*), intent(in)  :: filename
! Local Variables
integer :: i,j,n,u
! write vector to file
n = size(A,1)
open(newunit=u,file=filename,status="replace")
do i=1,n
    write(u,"("//FMT//")") A(i)
enddo
close(u)

end subroutine write_rvector