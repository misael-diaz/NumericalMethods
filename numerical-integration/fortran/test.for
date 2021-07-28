!
!   Applied Numerical Analysis                                July 27, 2021
!   ME 2020 FA21
!   Prof. M Diaz-Maldonado
!
!   source: test.for
!
!   Synopsis:
!   Tests the implementation of numerical integrators.
!
!
!   Copyright (C) 2021 Misael Diaz-Maldonado
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.


module objfun
    ! encapsulates the nonlinear function f(x) in a module for convenience
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    public
    contains
        function fun(x) result(f)
            ! Synopsis: Defines the nonlinear function, f(x)
            real(kind = real64), intent(in) :: x
            real(kind = real64):: f
            f = dexp(x)
            return
        end function
end module


program tests
    ! Synopsis:
    ! Integrates f(x) = exp(x) in the interval [a, b] using n intervals
    ! via Rectangle and Trapezoid Methods.
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use num_integrators, only: lsum, rsum, trap
    use objfun, only: fun
    implicit none
    procedure(fun), pointer :: fp => null()
    real(kind = real64) :: a, b, ni(3)
    integer(kind = int32), parameter :: n = 255

    fp => fun
    a = 0.0_real64
    b = 1.0_real64
    ni(1) = lsum (a, b, n, fp)        ! Left Riemann
    ni(2) = rsum (a, b, n, fp)        ! Right Riemann
    ni(3) = trap (a, b, n, fp)        ! Trapezoid Method

    print '(1X, A, F8.6)', "Left  Riemann:    ", ni(1)
    print '(1X, A, F8.6)', "Right Riemann:    ", ni(2)
    print '(1X, A, F8.6)', "Trapezoid Method: ", ni(3)
end program


! TODO:
! [ ] tabulate results as in the other tests.
