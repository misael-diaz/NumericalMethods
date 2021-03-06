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
            ! Synopsis:
            ! Defines the nonlinear function, f(x), to be integrated.
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
    real(kind = real64) :: a, b, ei, ni(3), err(3)
    integer(kind = int32), parameter :: n = 255

    fp => fun                   ! associates f[unction] p[ointer] to f(x)
    a = 0.0_real64              ! lower integration limit
    b = 1.0_real64              ! upper integration limit
    ei = fp(b) - fp(a)          ! exact integral for f(x)

    ni(1) = lsum (a, b, n, fp)  ! Left Riemann
    ni(2) = rsum (a, b, n, fp)  ! Right Riemann
    ni(3) = trap (a, b, n, fp)  ! Trapezoid Method

    err = abs(ei - ni) / ei * 100.0_real64

    print *, new_line('N')
    print '(A, 7X, A, 10X, A)', "Numerical Method", "Result", "% Error"
    print '(A)', "----------------------------------------------------"
    print '(A, F12.6, 8X, E9.3)', "Left  Riemann:    ", ni(1), err(1)
    print '(A, F12.6, 8X, E9.3)', "Right Riemann:    ", ni(2), err(2)
    print '(A, F12.6, 8X, E9.3)', "Trapezoid Method: ", ni(3), err(3)
    print *, new_line('N')
end program


! TODO:
! [x] tabulate results as in the other tests.
