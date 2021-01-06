!Module for the Differential evolution algorithm
module DE_mod
  implicit none
  private
  public :: DE_Method
contains
  !Differential evolution algorithm
  !input parameters
  !fob: objective function
  !m: population size
  !gmax: number of generations
  !bL, bU: limits for the search space (low and upper vertices of a
  !hyperrectangle)
  !Cr: crossover probability
  !F: mutation scale factor (ussually taken as 1)
  !op: variant of the method; op='D' is for the Dither variant, op='J'
  !is for the jitter variant (in our experience, the Dither variant has
  !better convergence)
  !
  !output parameters
  !minimo: the minimum
  !mejor: the best individual
  !tiempo: running time
  subroutine DE_Method(fob,m,gmax,bL,bU,Cr,F,op,minimo, mejor, tiempo)
    interface
      function fob(r) result(fobr)
        real(8), intent(in), dimension(:) :: r
        real(8) :: fobr
      end function fob
    end interface
  !
    integer, intent(in) :: m, gmax
    real(8), intent(in),dimension(:) :: bL, bU
    real(8), intent(in) :: Cr, F
    character(len=1), intent(in) :: op
    real(8), intent(out) :: minimo, tiempo
    real(8), intent(out), dimension(size(bL)) :: mejor
    real(8), dimension(m,size(bL)) :: X, V, U
    integer :: g
    real(8) :: t1, t2
    !
    call cpu_time(t1)
    X = X0(m, bL, bU)
    !
    do g = 1, gmax
      V = Vmut(X, bL, bU, F, op)
      U = Ucr(X, V, Cr)
      X = Sel(fob, X, U)
    end do
    !
    mejor = xbest(fob,X)
    minimo = fob(mejor)
    call cpu_time(t2)
    tiempo = t2 - t1
  end subroutine DE_Method
  !
  !function to extract from X the best individual according to the
  !objective function "fob"
  function xbest(fob, X) result(fr)
    interface
      function fob(r) result(fobr)
        real(8), intent(in), dimension(:) :: r
        real(8) :: fobr
      end function fob
    end interface
    !
    real(8), intent(in), dimension(:, :) :: X
    real(8), dimension(size(X, 2)) :: fr
    integer :: i, m
    real(8), dimension(size(X, 1)) :: fvals
    integer, dimension(1) :: minposv
    m = size(X, 1)
    do i=1, m
       fvals(i) = fob(X(i, :))
    end do
    minposv = minloc(fvals)
    fr = X(minposv(1), :)
  end function xbest
  !generates a random permutation from 1,2,...m of size nrs
  function randperm(m, nrs) result(fr)
    integer, intent(in) :: m, nrs
    integer, dimension(nrs) :: fr
    integer :: k, rnd
    integer, dimension(m) :: vec
    do k=1, m
       vec(k) = k
    end do
    do k=1, nrs
      rnd = floor(rand01() * (m - k + 1)) + 1
      fr(k) = vec(rnd)
      vec(rnd) = vec(m - k + 1)
    end do
  end function randperm
  !Selection operator
  function Sel(fob, X, U) result(Smat)
    interface
      function fob(r) result(fobr)
        real(8), intent(in), dimension(:) :: r
        real(8) :: fobr
      end function fob
    end interface
    real(8), intent(in), dimension(:,:) :: X, U
    real(8), dimension(size(X, 1), size(X, 2)) :: Smat
    integer :: i, m
    real(8), dimension(size(X,2)) :: ui, xi
    Smat = X
    m = size(X, 1)
    do i = 1, m
      xi = X(i, :)
      ui = U(i, :)
      if(fob(ui) < fob(xi)) Smat(i, :) = ui
    end do
  end function Sel
  !Crossover operator
  function Ucr(X, V, Cr) result(fr)
    real(8), intent(in), dimension(:,:) :: X, V
    real(8), intent(in) :: Cr
    real(8), dimension(size(x, 1), size(x, 2)) :: fr
    integer :: i, j, jrand, m, n
    real(8) :: randj
    m = size(X, 1)
    n = size(X, 2)
    fr = X
    do i=1, m
       do j=1, n
          randj = rand01()
          jrand = floor(rand01() * n) + 1
          if(randj <= Cr .or. j == jrand) fr(i, j) = V(i, j)
       end do
    end do
  end function Ucr
  !
  !mutation operator
  !D(J) is for the dither(jitter) variant
  function Vmut(X, bL, bU, F, op) result(fr)
    real(8), intent(in), dimension(:, :) :: X
    real(8), intent(in), dimension(size(X, 2)) :: bL, bU
    real(8), intent(in) :: F
    character(len=1), intent(in) :: op
    real(8), dimension(size(X,1),size(X,2)) :: fr
    integer :: i,j,m,n,r0,r1,r2
    integer, dimension(3) :: r
    real(8) :: rnd
    m = size(X,1)
    n = size(X,2)
    if(op == "D" .or. op == "d") then
       do i=1,m
          r = randperm(m,3)
          r0 = r(1)
          r1 = r(2)
          r2 = r(3)
          rnd = rand01()
          do j=1, n
             fr(i, j)  = X(r0,j) + F * rnd * (X(r1, j) - X(r2 ,j))
             if(fr(i, j) > bU(j) .or. fr(i, j) < bL(j)) fr(i, j) = bL(j) + &
                  rand01() * (bU(j) - bL(j))
          end do
       end do
    elseif(op == "J" .or. op == "j") then
       do i=1,m
          r = randperm(m,3)
          r0 = r(1)
          r1 = r(2)
          r2 = r(3)
          do j=1,n
             rnd = rand01()
             fr(i, j) = X(r0, j) + F * rnd * (X(r1, j) - X(r2, j))
             if(fr(i, j) > bU(j) .or. fr(i, j) < bL(j)) fr(i, j) = bL(j) + &
                  rand01() * (bU(j) - bL(j))
          end do
       end do
    else
       print*,"variant must be J or D"
    endif
  end function Vmut

  !Generates an initial population
  function X0(m, bL, bU) result(fr)
    integer, intent(in) :: m
    real(8), intent(in), dimension(:) :: bL, bU
    real(8), dimension(m, size(bL)) :: fr
    integer :: i,  j
    !
    do i = 1, m
       do j = 1, size(bL)
          fr(i, j) = bL(j) + rand01() * (bU(j) - bL(j))
       end do
    end do
  end function X0
  !generates a random number 0 y 1
  function rand01() result(fr)
    real(8) :: fr
    call random_number(fr)
  end function rand01
end module DE_mod
