!Runge-Kutta of 4th order
module RK4_mod
  implicit none

  private

  public :: RK4

contains
  
  !Runge-Kutta algorithm
  !Fvec is a vector function whose components correpond to the right
  !hand side of a first order differential equations system
  !t0, tf: limits of the domain
  !x0: initial conditions
  !npts: number of used intervals in the discretization of the domain
  !see below for the dimensions and kind of variables
  function RK4(Fvec,t0,tf,x0,npts) result(matf_r)
    interface
       function Fvec(t,x) result(Fvec_r)
         real(8), intent(in) :: t
         real(8), intent(in), dimension(:) :: x
         real(8), dimension(size(x)) :: Fvec_r
       end function Fvec
    end interface

    real(8), intent(in) :: t0, tf
    real(8), intent(in), dimension(:) :: x0
    integer, intent(in) :: npts

    !we included the initial point
    real(8), dimension(npts+1,size(x0)+1) :: matf_r
    real(8) :: h, ti, timas1
    real(8), dimension(npts+1) :: tvec
    real(8), dimension(npts+1,size(x0)) :: Xmat
    real(8), dimension(size(x0)) :: k1, k2, k3, k4
    real(8), dimension(size(x0)) :: xi, ximas1
    integer :: i, ncols

    h = (tf-t0)/npts
    tvec(1) = t0
    Xmat(1,:) = X0

    ncols = size(x0) + 1

    do i = 1, npts
       ti = tvec(i)
       xi = Xmat(i,:)

       k1 = h*Fvec(ti,xi)
       k2 = h*Fvec(ti + 0.5d0*h, xi + 0.5d0*k1)
       k3 = h*Fvec(ti + 0.5d0*h, xi + 0.5d0*k2)
       k4 = h*Fvec(ti + h, xi + k3)

       timas1 = ti + h
       ximas1 = xi + (k1 + 2d0*k2 + 2d0*k3 + k4)/6d0

       Xmat(i+1,:) = ximas1
       tvec(i+1) = timas1
    end do

    matf_r(:,1) = tvec
    matf_r(:,2:ncols) = Xmat 
  end function RK4
end module RK4_mod



