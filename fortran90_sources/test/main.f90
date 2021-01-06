include './modules/RK4_mod.f90'
include './modules/DE_mod.f90'
include './modules/tools_mod.f90'
! Load model parameters
module par_mod
    implicit none
    !------
    !{qr,eps}= {{0.6, 0.6},{0.6, 0.5},{0.6, 0.3},{0.5, 0.6},{0.5, 0.5},&
    !           {0.5, 0.3}, {0.3, 0.6}, {0.3, 0.5}, {0.3, 0.3}}
    !------
    real(8), parameter :: betaS = 0.36328215d0
    real(8), parameter :: betaA = 0.25152105d0
    real(8), parameter :: epsilon_par = 0.95d0
    real(8), parameter :: varepsilon = 0.1d0
    real(8), parameter :: kappa = 0.19607843d0
    real(8), parameter :: delta_h = 0.2d0
    real(8), parameter :: delta_v = 1d0/360d0
    real(8), parameter :: delta_R = 1d0/180d0
    real(8), parameter :: p = 0.1213d0
    real(8), parameter :: th = 0.2d0
    real(8), parameter :: gamma_a = 0.16750419d0
    real(8), parameter :: gamma_s = 0.09250694d0
    real(8), parameter :: gamma_h = 0.119869d0
    real(8), parameter :: mu = 3.913894d-5
    real(8), parameter :: mu_i_s = 0.0004d0
    real(8), parameter :: mu_h = 0.01632d0
    real(8), parameter :: lambda_l = 0.002777777777777778d0
    real(8), parameter :: lambda_v = 0.0006113521953813965d0
    real(8), parameter :: aD = 7.25d0
    real(8), parameter :: aIS = 0.0020127755438256486d0
    real(8), parameter :: aIH = 0.001411888738103725d0
    real(8), parameter :: c_L = 0.1
    real(8), parameter :: c_V = 0.01
    real(8), parameter :: xcoverage = 0.75d0
    real(8), parameter :: B = 0.000359216658d0, N = 26.446435d6
    real(8), parameter :: qr = 0.3d0
    !-------
    real(8), parameter :: L0 = 2.66260097d-1
    real(8), parameter :: S0 = 4.63606046d-1
    real(8), parameter :: E0 = 6.7033000d-4
    real(8), parameter :: IS0 = 9.28300000d-5
    real(8), parameter :: IA0 = 1.20986000d-3
    real(8), parameter :: H0 = 1.34157969d-4
    real(8), parameter :: R0 = 2.66125939d-1
    real(8), parameter :: D0 = 1.90074000d-3
    real(8), parameter :: V0 = 0.0d0
    real(8), parameter :: Xvac0 = 0.0d0
    real(8), parameter :: YIS0 = 0.12258164d0
    real(8), parameter :: YIH0 = 0.024516328000000004d0
    real(8), parameter :: XMAY0 = 0.0d0
    !-------
    real(8), parameter :: T0 = 0d0
    real(8), parameter :: Tf = 360d0
    integer, parameter :: ncon = 2
    integer, parameter :: npRK = Tf * 2, dimf = 13, npts = npRK + 1
    !
    real(8), parameter, dimension(dimf) :: x0 =&
        [L0, S0, E0, IS0, IA0, H0, R0, D0, V0, YIS0, YIH0, Xvac0, XMAY0]
 end module par_mod
  !
  !objective function module
module fob_mod
    use par_mod
    use tools_mod
    use RK4_mod
    implicit none
    private
    public :: fob, qj
    contains
    !since we takie controls of equal dimensions, x must be a
    !multiple of ncon
        function fob(x) result(fob_r)
            real(8), intent(in), dimension(:) :: x
            real(8) :: fob_r
            real(8), dimension(size(x) / ncon) :: rL, rV
            real(8), dimension(npts) :: ct, cL, cS, cE, cIs, cIa, cH, cD, cR
            real(8), dimension(npts) :: cYIS, cYIH, cXvac, CXMAY, cuL, cuV
        !    real(8), dimension(npts) :: integ, integ_aux
        !    real(8), dimension(npts, 2) :: matinteg, matinteg_aux
            real(8), dimension(npts, dimf + 1) :: Amat
            integer :: i, dimx, dimr
            real(8) :: DT, xvacT, XMAYT, YIST, YIHT
        !
            dimx = size(x)
            dimr = dimx / ncon
            rL = x(1:dimr)
            rV = x(dimr + 1:ncon * dimr)
        !
             Amat = RK4(Fvec, T0, Tf, x0, npRK)
        !   x0 = [L0,S0,E0,IS0,IA0,H0,R0,D0,V0] are the initial conditions
        !   Extracting the columns for the integral. The first column
        !   corresponds to time.
            ct  = Amat(:, 1)
            cL  = Amat(:, 1 + 1)
            cS  = Amat(:, 2 + 1)
            cE  = Amat(:, 3 + 1)
            cIs = Amat(:, 4 + 1)
            cIa = Amat(:, 5 + 1)
            cH = Amat(:, 6 + 1)
            cR  = Amat(:, 7 + 1)
            cD =  Amat(:, 8 + 1)
            cYIS = Amat(:, 10+1)
            cYIH = Amat(:, 11+1)
            cXvac = Amat(:, 12+1)
            cXMAY = Amat(:, 13 + 1)
        !
            fob_r = 1d10
            if(kappa * maxval(cIs) > B) return
            !
            do i =1, npts
                cuV(i) = qj(ct(i), rV)
                cuL(i) = qj(ct(i), rL)
            end do
            ! Vacine doses counter
            ! integ_aux = (cuV + lambda_v) * (cL + cS + cE + cIa + cR)
            ! matinteg_aux(:, 1) = ct
            ! matinteg_aux(:, 2) = integ_aux
            !
            !xvacT = trapz(matinteg_aux)
            xvacT = cXvac(ubound(cXvac, dim=1))
            if(xvacT < xcoverage) return
            if(xvacT > .90d0) return
            !integ = aIS * p * kappa * cE + &
            ! (aD * mu_i_s + aH * delta_h) * cIs + &
            !    cuL ** 2 + cuV ** 2
            !
            !The integrand as a set of points (t,integrand(t))
            !matinteg(:, 1) = ct
            !matinteg(:, 2) = integ
            !fob_r = trapz(matinteg)
            !
            XMAYT = cXMAY(ubound(cXMAY, dim=1))
            YIST = cYIS(ubound(cYIS, dim=1))
            DT = cD(ubound(cYIS, dim=1))
            fob_r = XMAYT + aIS * (YIST - YIS0) +  aD * (DT - D0) + &
                    aIH * (YIHT - YIH0)
            contains
            !The Fvec function for the SED
                function Fvec(t, xvec) result(vecf_r)
                real(8), intent(in) :: t
                real(8), intent(in), dimension(:) :: xvec
                real(8), dimension(size(xvec)) :: vecf_r
                real(8) :: xL, xS, xE, xIS, xIA, xH, xR, xD, xV
                real(8) :: yIS, yIH, Vcon, XMAY
                real(8) :: fL, fS, fE, fIS, fIA, fH
                real(8) :: fR, fD, fV, fYIS, fYIH, fVcon, fXMAY, u_Lt, u_Vt
                real(8) :: Nstar, lambda
                !
                !x = [L,S,E,IS,IA,H,R,D,V]
                xL  = xvec(1)
                xS  = xvec(2)
                xE  = xvec(3)
                xIS = xvec(4)
                xIA = xvec(5)
                xH  = xvec(6)
                xR  = xvec(7)
                xD  = xvec(8)
                xV  = xvec(9)
                yIS = xvec(10)
                yIH = xvec(11)
                Vcon = xvec(12)
                XMAY = xvec(13)
                !
                Nstar = xL + xS + xE + xIS + xIA + xH + xR + xV
                !
                lambda = (betaS * xIS + betaA * xIA ) / Nstar
                u_Lt = qj(t, rL)
                u_Vt = qj(t, rV)
                fL = th * mu * Nstar - &
                    ( varepsilon * lambda + &
                     lambda_l + u_Lt + &
                     lambda_v + u_Vt + &                    
                      mu) * xL
                fS = (1d0 - th) * mu * Nstar + & 
                    (lambda_l + u_Lt) * xL + &
                    delta_v * xV + delta_R * xR - &
                    (lambda + lambda_v + u_Vt + mu) * xS
                fE = lambda * &
                        (varepsilon * xL + (1d0 - epsilon_par) * xV + xS) - &
                        (kappa + mu) * xE
                fIS = p * kappa * xE - (gamma_s + mu_i_s + delta_h + mu) * xIS
                fIA = (1d0 - p) * kappa * xE - (gamma_a + mu) * xIA
                fH = delta_h * xIS - (gamma_h + mu_h + mu) * xH
                fR = gamma_s * xIS + gamma_a * xIA + gamma_h*xH - &
                    (delta_R + mu)*xR
                fD = mu_i_s * xIS + mu_h * xH
                fV = (lambda_v +  u_Vt) * (xS + xL) - &
                    ((1d0 - epsilon_par) * lambda +  delta_v + mu) * xV
                fYIS = p * kappa * xE
                fYIH = delta_h * xIS
                fVcon =  (lambda_v +  u_Vt) * (xL + xS + xE + xIA +xR)
                fXMAY = 0.5 * c_L * u_Lt ** 2 + 0.5 * c_V * u_Vt **2
                !the right hand side of the SED
                vecf_r = [fL, fS, fE, fIS, fIA, fH, fR, fD, fV, &
                            fYIS, fYIH, fVcon, fXMAY]
            end function Fvec
        end function fob
    !
    !the controls
    !auxiliar function for fob:  qj, equally spaced intervals
    !modify this functions if the intervals are of different legths
    !for more informatio see:
    !Engineering Applications of Artificial Intelligence 97 (2021)
    !104086; https://doi.org/10.1016/j.engappai.2020.104086
        function qj(t, r) result(fr)
            real(8), intent(in) :: t
            real(8), intent(in), dimension(:) :: r
            real(8) :: fr
            real(8) :: dt, eps
            integer :: idt
            eps = epsilon(eps)
            dt = (Tf - T0) / size(r)
            idt = floor(t / dt) + 1
            if(abs(t - Tf) <= eps) idt = size(r)
            fr = r(idt)
        end function qj
end module fob_mod
!
!main program; the optimization
!
program main
    use RK4_mod
    use tools_mod
    use fob_mod
    use DE_mod
    use par_mod, only : T0, Tf, npRK, ncon,  lambda_v, lambda_l
    implicit none
    !
    !dimx must be multiple of  ncon
    integer, parameter :: dimx = Tf / 2, dimr = dimx / ncon, nt = npRK + 1
    integer, parameter :: m = dimx * 4, gmax = 3000
    real(8), parameter :: F = 1d0, Cr = 0.3d0
    character(len = 1), parameter :: JD = "D"
    !real(8), parameter, dimension(dimx) :: bL = -.05d0 * lambda_v
    !real(8), parameter, dimension(dimx) :: bU = 5.0d0 * lambda_v
    real(8), dimension(dimx) :: bL
    real(8), dimension(dimx) :: bU
    real(8), dimension(dimx) :: mejor
    real(8) :: minimo, tiempo, ti
    real(8), dimension(dimr) :: rL, rV
    character(len=:), allocatable :: file_con
    real(8), dimension(nt) :: tvec
    real(8), dimension(nt, ncon + 1) :: controles !the first column is time
    integer :: i
    !
    !Control bounds
    ! Lower bound for Lockdown and vaccination controls
    bL(1:dimr) = -4.0 * lambda_L 
    bL(dimr+1:dimx) = -lambda_v
    !Upper boudn of Lockdown and vaccination controls
    bU(1:dimr) = 3d0 * lambda_L
    bU(dimr+1:dimx) = 7d0 * lambda_v
    call DE_Method(fob ,m, gmax, bL, bU, Cr, F, JD, minimo, mejor, tiempo)
    !
    !The best individual
    open(unit = 11, file = 'best.dat', status = 'unknown')
    do i=1, size(mejor)
       write(11,"(*(f0.8, 2x))") mejor(i)
    end do
    close(unit = 11)
    !-----------
    rL = mejor(1:dimr)
    rV = mejor(dimr + 1:ncon * dimr)
    !The controls
    file_con = 'rvecs_'//IntToString(dimx)//'_'//IntToString(m)// &
         '_'//IntToString(gmax)//'.dat'
    open(unit = 11, file = file_con, status = "unknown")
    do i=1, size(rL)
       write(11, "(*(f0.8,2x))") rL(i), rV(i)
    end do
    close(unit = 11)
    !-----------
    tvec = linspace(T0, Tf, nt)
    do i=1, nt
       ti = tvec(i)
       controles(i, :) = [ti, qj(ti, rL), qj(ti, rV)]
    end do
    ![t,uL(t),uV(t)] the controls evaluated in tvec
    file_con = 'datos_'//IntToString(dimx)//'_'//IntToString(m)//      &
         '_'//IntToString(gmax)//'.dat'
    open(unit = 11, file = file_con, status = "unknown")
    do i = 1, nt
       write(11, "(*(f0.8,2x))") controles(i, :)
    end do
    close(unit = 11)
    !-----------
    print*, 't(s)=',tiempo
    print*,'min=',minimo
end program
