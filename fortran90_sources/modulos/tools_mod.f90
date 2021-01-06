module tools_mod
  implicit none
  private

  public :: trapz, IntToString, linspace

contains

  !trapz function
  !Trapezoidal numerical integration
  !X: matrix of points [x,f(x)]
  function trapz(X) result(fr)
    real(8), intent(in), dimension(:,:) :: X
    real(8) :: fr
    real(8), dimension(size(X,1)) :: xv, yv
    integer :: n

    n = size(X,1)

    xv = X(:,1)
    yv = X(:,2)

    fr =  0.5d0*sum((xv(2:n) - xv(1:n - 1))*(yv(1:n - 1) + yv(2:n)))
  end function trapz

  !auxiliar function for changing an integer to string
  function IntToString(n) result(f_result)
    character(len = :), allocatable :: f_result
    integer, intent(in) :: n
    character(len = 100) :: auxf_r

    write(auxf_r,*) n
    f_result = trim(adjustl(auxf_r))
  end function IntToString

  !Generate linearly spaced vector
  !a: initial point
  !b: final point
  !n: mumber of points
  function linspace(a,b,n) result(fr)
    real(8), intent(in) :: a, b
    integer, intent(in) :: n
    real(8), dimension(n) :: fr
    real(8) :: dx
    integer :: k

    if(n==1) then
       fr(1) = a
       return
    end if

    dx = (b-a)/(n-1)    

    do k=1, n
       fr(k) = a + (k-1)*dx
    end do
  end function linspace

end module tools_mod
