  ! to compile:
  ! GNU/Linux
  ! ifort -heap-arrays 1024 InfectedCovid.f90
  !
  ! or using the gfortran compiler:
  ! gfortran -O3 optimiza.f90
  !
  ! Windows
  ! ifort /F:2147483647
  !
  ! to execute
  ! GNU/Linux
  ! ./*.out
  !
  ! Windows
  ! *.exe
  !
  ! '*' stands for the name of the executable created after compilation
  
  !-------
  !npRK is the number of intervals used in the discretization for the
  !RK4 method
  !npts = npRK+1 since the initial point is included
  
  include '../modulos/tools_mod.f90'
  include '../modulos/RK4_mod.f90'

  !module for parameters
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
    real(8), parameter :: delta_v = 1d0/120d0 
    real(8), parameter :: delta_R = 1d0/180d0
    real(8), parameter :: p = 0.1213d0
    real(8), parameter :: th = 0.2d0
    real(8), parameter :: gamma_a = 0.16750419d0
    real(8), parameter :: gamma_s = 0.09250694d0
    real(8), parameter :: gamma_h = 5.079869d-1
    real(8), parameter :: mu = 3.913894d-5    
    real(8), parameter :: mu_i_s = 0.01632d0
    real(8), parameter :: mu_h = 0.01632d0
    real(8), parameter :: lambda_v = 0.0006113521953813965d0
    real(8), parameter :: aD = 1d0
    real(8), parameter :: aIS = 1d0, aH = 1d0    
    real(8), parameter :: xcoverage = 0.2d0
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
    !-------
    real(8), parameter :: T0 = 0d0
    real(8), parameter :: Tf = 360d0
    integer, parameter :: ncon = 2
    integer, parameter :: npRK = Tf*1, dimf = 9, npts = npRK+1

    real(8), parameter, dimension(dimf) :: x0 = [L0,S0,E0,IS0,IA0,H0,  &
         R0,D0,V0]

  end module par_mod

  !modulo para ClasesMat y qj
  module dynsol_mod
    use par_mod
    use tools_mod
    use RK4_mod
    implicit none
    private
    public :: ClasesMat, qj

  contains

    !la dimension de x debe ser un multiplo de 4, tomaremos los
    !controles de las mismas dimensiones
    function ClasesMat(x) result(Amat)
      real(8), intent(in), dimension(:) :: x
      real(8), dimension(npts, dimf + 1) :: Amat 
      real(8), dimension(size(x)/ncon) :: rL, rV
      integer :: dimx, dimr
      real(8) :: xvacT

      dimx = size(x)
      dimr = dimx/ncon
      rL = x(1:dimr)
      rV = x(dimr+1:ncon*dimr)

      Amat = RK4(Fvec,T0,Tf,x0,npRK)

    contains
      !The Fvec function for the SED
      function Fvec(t, xvec) result(vecf_r)
        real(8), intent(in) :: t
        real(8), intent(in), dimension(:) :: xvec
        real(8), dimension(size(xvec)) :: vecf_r
        real(8) :: xL, xS, xE, xIS, xIA, xH, xR, xD, xV
        real(8) :: fL, fS, fE, fIS, fIA, fH, fR, fD, fV
        real(8) :: Nstar, lambda

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

        Nstar = xL + xS + xE + xIS + xIA + xH + xR + xV

        lambda = (betaS*xIS + betaA*xIA)/Nstar

        fL = th*mu*Nstar - epsilon_par*lambda*xL - qj(t,rL)*xL - mu*xL

        fS = (1d0 - th)*mu*Nstar + qj(t,rL)*xL + delta_v*xV + delta_R* &
             xR - (lambda + lambda_v + qj(t,rV) + mu)*xS

        fE = lambda*(epsilon_par*xL + (1d0 - varepsilon)*xV + xS) -    &
             (kappa + mu)*xE

        fIS = p*kappa*xE - (gamma_s + mu_i_s + delta_h + mu)*xIS

        fIA = (1d0 - p)*kappa*xE - (gamma_a + mu)*xIA

        fH = delta_h*xIS - (gamma_h + mu_h + mu)*xH 

        fR = gamma_s*xIS + gamma_a*xIA + gamma_h*xH - (delta_R + mu)*xR

        fD = mu_i_s*xIS + mu_h*xH 

        fV = (lambda_v +  qj(t,rV))*xS - ((1d0 - varepsilon)*lambda +  &
             delta_v + mu)*xV

        !the right hand side of the SED
        vecf_r = [fL, fS, fE, fIS, fIA, fH, fR, fD, fV]
      end function Fvec
    end function ClasesMat

    !los controles
    !funcion auxiliar para fob:  qj, equally spaced intervals
    !modify this functions if the intervals are of different legths
    function qj(t, r) result(fr)
      real(8), intent(in) :: t
      real(8), intent(in), dimension(:) :: r
      real(8) :: fr
      real(8) :: dt, eps
      integer :: idt

      eps = epsilon(eps)
      dt = (Tf - T0)/size(r)

      idt = floor(t/dt) + 1
      if(abs(t-Tf) <= eps) idt = size(r)

      fr = r(idt)
    end function qj
  end module dynsol_mod

  !main program
  program test
    use RK4_mod 
    use tools_mod
    use dynsol_mod
    use par_mod, only : npts, dimf, Tf
    implicit none

    !ver dimx en el archivo optimiza.f90 debe ser la misma dimx de aqui.
    !mejor es generado al ejecutar optimiza.f90
    integer, parameter :: dimx = Tf/2 
    real(8), dimension(dimx) :: mejor, cero
    real(8), dimension(npts, dimf + 1) :: AmatCV
    integer :: i, unit_num
    character(len=:),dimension(:), allocatable :: header
    integer :: decimal_places

    call unit_assign(unit_num)

    open(unit = unit_num, file = 'mejor.dat', status = 'old')
    do i = 1, dimx
       read(unit_num,*) mejor(i)
    end do
    close(unit = unit_num)

    !controlled case
    AmatCV = ClasesMat(mejor)
    header = ['CL','CS','CE','CIS','CIA','CH','CR','CD','CV']
    decimal_places = 9
    call savemat2csv(AmatCV,'ControlClases.csv',header,decimal_places)

    !No controlled case
    cero = 0.0d0
    AmatCV = ClasesMat(cero)
    header = ['NCL','NCS','NCE','NCIS','NCIA','NCH','NCR','NCD','NCV']
    call savemat2csv(AmatCV,'NoControlClases.csv',header,decimal_places)
  end program test

