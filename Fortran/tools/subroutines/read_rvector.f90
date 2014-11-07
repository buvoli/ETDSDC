! Reads size(A) lines from specificed file, written by write_rvector function

subroutine read_rvector(A,filename)

! Arguments
character(len=*), intent(in) :: filename
real(dp),intent(out) :: A(:)
! Local Variables
real(dp) :: val
integer :: i,u,nA
character(len=*), parameter  :: FMT = "ES23.15E3"
! Read Vector
nA = size(A)
open(newunit=u,file=filename, status="old")
do i=1,nA
    read(u,"("//FMT//")") val
    A(i) = val
enddo
close(u)

end subroutine read_rvector