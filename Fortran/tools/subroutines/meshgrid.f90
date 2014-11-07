subroutine meshgrid(xr,yr,x,y)
! Replicates the grid vectors xr and yr to produce a full grid represented by the output coordinate arrays X and Y.
! === Parameters ===
! xr - (vector) x grid points
! yr - (vector) y grid points
! === Output ===
! x  - (matrix) x values at gridpoints
! y  - (matrix) y values at gridpoint
real(kind=8), intent(in)  :: xr(:),yr(:)
real(kind=8), intent(out) :: x(size(xr),size(yr)),y(size(xr),size(yr))
integer i,j,n,m
n = size(xr)
m = size(yr)
do i=1,n
    do j=1,n
        x(i,j) = xr(j)
        y(i,j) = yr(i)
    enddo
enddo
end subroutine meshgrid