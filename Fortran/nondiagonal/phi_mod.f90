! ============================ Module Description ===========================
! PHI_MOD: Initilizes Phi Functions for matrix arguments
!   1. phi - initilizes phi functions using taylor-series and scaling and 
!            squaring integral & recursion.
! ===========================================================================

module phi_mod

    ! Module Parameters
    use tools_mod, only: dp, norm
    implicit none

    contains

    ! =======================================================================
    ! PHI   Evaluates \phi_i(L) for i=0,...,n and matrix L using
    !       Taylor series and scaling and squaring relation.
    !
    ! Arguments
    !
    !   L   (input) COMPLEX*16 array
    !       KxK matrix cooresponding to PDE linear component
    !
    !   n   (input) INTEGER
    !       highest phi function to initialize.
    !
    !   P   (output) COMPLEX*16
    !       KxKx(n+1) array where on exit P(:,:,i) = \phi_{i-1}(L)
    ! =======================================================================

    subroutine phi(L,n,P)

        ! Arguments
        complex(dp), intent(in)     :: L(:,:)
        integer,     intent(in)     :: n
        complex(dp), intent(out)    :: P(size(L,1),size(L,2),n+1)

        ! Variable Declarations
        integer,        parameter   :: M = 20                       ! Number of Taylor Terms to use

        real(dp)                    :: fct(0:M+n)                   ! Array to store factorials
        complex(dp), allocatable    :: PT(:,:,:), T(:,:), TL(:,:)   ! temporary storage arrays
        complex(dp), allocatable    :: L0(:,:)                      ! temporary storage of L/2^s
        integer                     :: i,j,k,nL,s
        real(kind=8)                :: sf

        nL = size(L,1)
        allocate(PT(nL,nL,n+1), T(nL,nL))

        ! Choosing scaling factor
        s  = max(0,ceiling(log(max(norm(L,1),norm(L,-1)))/log(2.d0)))
        L0 = L/2.0_dp**s

        ! Precompute Factorials
        fct(0) = 1.0_dp
        do i=1,M+n
            fct(i) = fct(i-1)*real(i,dp)
        enddo

        ! Initial Taylor Series Definition using Horner Scheme
        do i=0,n
            P(:,:,i+1) = L0/fct(M+i)
            do j=(M-1),1,-1
                ! Using MATMUL (slower)
                !P(:,:,i+1) = L0/fct(j+i) + matmul(L0, P(:,:,i+1)) !(USING MATMUL)
                ! Using BLAS
                TL = L0
                call ZGEMM('N','N',nL,nL,nL,COMPLEX(1.0_dp,0.0_dp),P(:,:,i+1),nL,L0,nL,COMPLEX(1.0_dp/fct(j+i), 0.0_dp),TL,nL)
                P(:,:,i+1) = TL
            enddo
            do j=1,nL
                P(j,j,i+1) = 1.0_dp/fct(i) + P(j,j,i+1)
            enddo
        enddo

        ! Scale Phi function by powers of two
        do i=1,s
            PT = P
            do j=0,n
                if(j == 0) then ! square exponential
                    ! Using MATMUL (slower)
                    !PT(:,:,1) = matmul(PT(:,:,1),PT(:,:,1))
                    ! Using BLAS
                    call ZGEMM('N','N',nL,nL,nL,COMPLEX(1.0_dp,0.0_dp),PT(:,:,1),nL,PT(:,:,1),nL,COMPLEX(0.0_dp,0.0_dp),T,nL)
                    PT(:,:,1) = T
                else  ! apply recrusion to higher phi functions
                    sf = 2.0_dp**(-j)
                    !PT(:,:,j+1) = sf * matmul(P(:,:,1),P(:,:,j+1))
                    call ZGEMM('N','N',nL,nL,nL,COMPLEX(sf,0.d0),P(:,:,1),nL,P(:,:,j+1),nL,COMPLEX(0.d0,0.d0),T,nL)
                    PT(:,:,j+1) = T
                    do k=1,j
                        PT(:,:,j+1) = PT(:,:,j+1) + (sf/fct(j-k)) * P(:,:,k+1)
                    enddo
                endif
            enddo
            P = PT
        enddo

        ! Deallocate Variables
        deallocate(PT,L0)

    end subroutine phi

end module phi_mod