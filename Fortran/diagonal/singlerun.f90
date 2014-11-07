! ============================ Program Description ==========================
! SINGLERUN runs ETDRK4, ETDSDC, IMEXSDC Methods on specified equation with 
!           fixed timestep for each method for PDEs with a nondiagonal linear
!           operators.
! ===========================================================================

program singlerun

    ! ===================== Specify Equation Module Here ====================
    use kuramoto_mod,     only: L,N,init,y0,Np,tspan ! kuramoto_mod, nikolaevskiy_mod, quasigeostrophic_mod, kdv_mod, nls_mod
    ! =======================================================================

    use tools_mod,   only: dp, chebpts, print_cvector
    use etdsdc_mod,  only: etdsdc, etdsdc_settings
    use imexsdc_mod, only: imexsdc, imexsdc_settings
    use etdrk4_mod,  only: etdrk4

    implicit none

    ! Local Variables
    complex(dp), allocatable    :: Lambda(:),N0(:)                                          ! PDE Linear Operator
    integer                     :: Nt                                                       ! Number of Timesteps
    complex(dp), allocatable    :: yr_etdsdc(:), yr_imexsdc(:), yr_etdrk(:), yr_etdmp(:)    ! Solution Storage Arrays
    real(dp)                    :: times(4)                                                 ! Time Arrays
    integer                     :: Npr                                                      ! Number of outputs to print

    ! Method Settings (See Method Modules for Type Descriptions)
    TYPE(etdsdc_settings)       :: etdsdc_s
    TYPE(imexsdc_settings)      :: imexsdc_s

    allocate(lambda(Np),N0(Np))
    allocate(yr_etdsdc(Np), yr_imexsdc(Np), yr_etdrk(Np))

    ! Initialize Equation Module & Set Linear Operator lambda
    call init()
    call L(lambda)
    Npr = min(10,size(y0)) ! Number of outputs to print

    ! == Useful for testing L and N functions ==
    !call print_cvector(y0(1:Npr),"y0")
    !call print_cvector(lambda(1:Npr),"L")
    !call N(0.0_dp,y0,N0)
    !call print_cvector(N0(1:Npr),"N0")

    ! === ETDSDC Sample Call ===
    Nt = 300
    etdsdc_s%tau = chebpts(8,[0.0_dp, 1.0_dp])
    etdsdc_s%m = 7
    call etdsdc(lambda,N,tspan,y0,Nt,etdsdc_s,yr_etdsdc,times)

    write (*,*) "=== ETDSDC Method ==="
    write (*,*) "CPU Time: ", times(1)
    write (*,*) "CLK Time: ", times(2)
    write (*,*) "=== ETDSDC Coefficients ==="
    write (*,*) "CPU Time: ", times(3)
    write (*,*) "CLK Time: ", times(4)
    call print_cvector(yr_etdsdc(1:Npr),"y_final")

    ! === IMEXSDC Sample Call ===
    Nt = 2000
    imexsdc_s%tau = chebpts(8,[0.d0, 1.d0])
    imexsdc_s%m = 7
    call imexsdc(lambda,N,tspan,y0,Nt,imexsdc_s,yr_imexsdc,times)

    write (*,*) "=== IMEXSDC Method ==="
    write (*,*) "CPU Time: ", times(1)
    write (*,*) "CLK Time: ", times(2)
    write (*,*) "=== IMEXSDC Coefficients ==="
    write (*,*) "CPU Time: ", times(3)
    write (*,*) "CLK Time: ", times(4)
    call print_cvector(yr_imexsdc(1:Npr),"y_final")

    ! === ETDRK4 Sample Call ===
    Nt = 20000
    call etdrk4(lambda,N,tspan,y0,Nt,yr_etdrk,times)

    write (*,*) "=== ETDRK4 Method ==="
    write (*,*) "CPU Time: ", times(1)
    write (*,*) "CLK Time: ", times(2)
    write (*,*) "=== ETDRK4 Coefficients ==="
    write (*,*) "CPU Time: ", times(3)
    write (*,*) "CLK Time: ", times(4)
    call print_cvector(yr_etdrk(1:Npr),"y_final")

end program singlerun