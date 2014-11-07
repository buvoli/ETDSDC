! ============================ Program Description ==========================
! SINGLERUN generates results for a numerical experiment using ETDRK4, ETDSDC,
!           IMEXSDC Methods on specified equation and saves data in results
!           folder.
! ===========================================================================
program experimentrun

    ! ===================== Specify Equation Module Here ====================
    use kdv_mod,     only: L,N,init,y0,Np,tspan,Fs,reference_methods,eqn_name, error_filter
    ! =======================================================================

    use tools_mod,   only: dp,chebpts, isfinite, write_rmatrix, write_rvector, write_cvector
    use etdrk4_mod,  only: etdrk4
    use etdsdc_mod,  only: etdsdc, etdsdc_settings
    use imexsdc_mod, only: imexsdc, imexsdc_settings

    implicit none

    ! Local Variables
    complex(kind=8),allocatable :: y_reference(:), y_out(:),Lambda(:,:)
    integer                     :: F,relcost,nn,Nt,solution_count,i,j,n_imex, n_etd, n_Fs, sdc_reference_n
    TYPE(etdsdc_settings)       :: etdsdc_s
    TYPE(imexsdc_settings)      :: imexsdc_s
    integer, allocatable        :: ETD_methods(:), IMEX_methods(:)
    real(kind=8)                :: time(4)
    real(kind=8),allocatable    :: times(:,:), errors(:,:), Nts(:,:), hs(:,:)
    character(len=*), parameter :: results_dir = "../results/"

    ! IMEXSDC and ETDSDC Methods to test
    ETD_methods  = [4, 8, 16, 32]
    IMEX_methods = [4, 8, 16, 32]

    ! Init Equation
    allocate(Lambda(Np,Np))
    call init()
    call L(Lambda)

    write (*,"(a,a,a)") "Running experiment for ", eqn_name, " equation"
    write (*,"(a,$)") "Computing reference solution... "

    ! === Compute Reference Solution ===
    allocate(y_out(Np),y_reference(Np))
    F = 4*Fs(size(Fs))
    y_reference = (0.d0,0.d0)
    solution_count = 0
    sdc_reference_n = 32

    if(reference_methods(1)) then ! Use ETDSDC_N^N-1
        etdsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
        etdsdc_s%m   = sdc_reference_n-1
        relcost      = size(etdsdc_s%tau) * (etdsdc_s%m + 1)
        Nt = F/relcost + 1
        call etdsdc(Lambda,N,tspan,y0,Nt,etdsdc_s,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    if(reference_methods(2)) then ! Use IMEXSDC_N^N-1
        imexsdc_s%tau = chebpts(sdc_reference_n,(/0.d0, 1.d0 /))
        imexsdc_s%m   = sdc_reference_n-1
        relcost       = size(imexsdc_s%tau) * (imexsdc_s%m + 1)
        Nt = F/relcost + 1
        call imexsdc(lambda,N,tspan,y0,Nt,imexsdc_s,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    if(reference_methods(3)) then ! Use ETDRK4
        relcost = 4
        Nt = F/relcost + 1
        call etdrk4(Lambda,N,tspan,y0,Nt,y_out,time)
        if(isfinite(y_out)) then
            y_reference = y_reference + y_out
            solution_count = solution_count + 1
        endif
    endif

    y_reference = y_reference/real(solution_count,8)
    write (*,"(a)") "done."

    ! === Run Numerical tests ===
    write (*,"(a)") "Running Numerical Tests... "
    n_etd  = size(ETD_methods)
    n_imex = size(IMEX_methods)
    n_Fs   = size(Fs)
    allocate(times(n_Fs,n_etd+n_imex+1),errors(n_Fs,n_etd+n_imex+1),Nts(n_Fs,n_etd+n_imex+1),hs(n_Fs,n_etd+n_imex+1))
    do i=1,n_Fs
        F = Fs(i)
        write (*,"(a,I4,a,I4,a,I10)") "     F (",i,"/",n_Fs,") = ",F
        ! Run ETD Methods
        do j=1,size(ETD_methods)
            nn = ETD_methods(j)
            etdsdc_s%tau = chebpts(nn,(/0.d0, 1.d0 /))
            etdsdc_s%m   = nn-1
            relcost       = nn**2
            Nt = F/(relcost) + 1
            call etdsdc(Lambda,N,tspan,y0,Nt,etdsdc_s,y_out,time)
            Nts(i,j)    = Nt
            errors(i,j) = error_filter(y_out,y_reference)
            times(i,j)  = time(2)
        enddo
        ! Run IMEX Methods
        do j=1,size(IMEX_methods)
            nn = IMEX_methods(j)
            imexsdc_s%tau = chebpts(nn,(/0.d0, 1.d0 /))
            imexsdc_s%m   = nn-1
            relcost       = nn**2
            Nt = F/(relcost) + 1
            call imexsdc(Lambda,N,tspan,y0,Nt,imexsdc_s,y_out,time)
            Nts(i,j+n_etd)    = Nt
            errors(i,j+n_etd) = error_filter(y_out,y_reference)
            times(i,j+n_etd)  = time(2)
        enddo
        ! Run RK4 Method
        relcost = 4
        Nt = F/(relcost) + 1
        call etdrk4(Lambda,N,tspan,y0,Nt,y_out,time)
        Nts(i,1+n_etd+n_imex)    = Nt
        errors(i,1+n_etd+n_imex) = error_filter(y_out,y_reference)
        times(i,1+n_etd+n_imex)  = time(2)
    enddo
    write (*,"(a)") "done."
    ! Store stepsizes
    hs = (tspan(2) - tspan(1))/Nts

    ! Save Results
    call write_rmatrix(Nts,     results_dir//eqn_name//"/Nts.txt")
    call write_rmatrix(hs,      results_dir//eqn_name//"/hs.txt")
    call write_rmatrix(errors,  results_dir//eqn_name//"/errors.txt")
    call write_rmatrix(times,   results_dir//eqn_name//"/times.txt")
    call write_rvector(Fs,      results_dir//eqn_name//"/Fs.txt")
    call write_cvector(y_reference,      results_dir//eqn_name//"/reference.txt")

end program experimentrun