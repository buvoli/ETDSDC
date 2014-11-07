function norm(A,nrm_type) result(nrm)
! Returns matrix norm of A. Only Supports 1 norm and Inf norm
! === Parameters ===
! A         - (matrix)
! nrm_type  - (integer) nrm_type=1 for one norm and nrm_type=-1 for Inf norm
! === Output ===
! nrm       - (double) matrix norm
complex(kind=8),  intent(in) :: A(:,:)
integer, intent(in)          :: nrm_type
real(kind=8)                 :: nrm, nrm_canidate
integer :: i,j,d1,d2
d1 = size(A,1)
d2 = size(A,2)
nrm = 0.d0
if(nrm_type == -1) then     ! Infinity nrm
    do i=1,d1
        nrm_canidate = 0.d0
        do j=1,d2
            nrm_canidate = nrm_canidate + abs(A(i,j))
        enddo
        if(nrm_canidate > nrm) then
            nrm = nrm_canidate
        endif
    enddo
elseif(nrm_type == 1) then ! 1 nrm
    do j=1,d2
        nrm_canidate = 0.d0
        do i=1,d1
            nrm_canidate = nrm_canidate + abs(A(i,j))
        enddo
        if(nrm_canidate > nrm) then
            nrm = nrm_canidate
        endif
    enddo
endif
end function norm