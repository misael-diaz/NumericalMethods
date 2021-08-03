!
!   Applied Numerical Analysis                                July 26, 2021
!   ME 2020 FA21
!   Prof. M Diaz-Maldonado
!
!   source: test_solvers.for
!
!   Synopsis:
!   Tests nonlinear solvers.
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

            !< Note the use of `real64' for values with double precision >!
            f = ( 1.0_real64 / sqrt(x) + 2.0_real64 *  &
                & log10(0.024651_real64 / 3.7_real64 + &
                & 2.51_real64 / (9655526.5_real64 * sqrt(x) ) ) )

            return
        end function
end module


program tests
    use, intrinsic :: iso_fortran_env, only: real64
    use nlsolvers, only: bisect, regfal, nls_conf
    use objfun, only: fun
    implicit none
    type(nls_conf) :: opts
    procedure(fun), pointer :: fp => null()
    real(kind = real64) :: lb, ub
    real(kind = real64) :: x1, x2

    ! sets the tolerance and maximum number of iterations of the solver
    opts % tol      = 1.0e-12_real64
    opts % max_iter = 256

    ! defines the (root) bracketing interval [lb, ub]
    lb = 1.0e-2_real64
    ub = 9.0e-2_real64


    fp => fun   !< associates procedure pointer to f(x) >!


    ! solves for the root of f(x) numerically
    x1 = bisect (lb, ub, fp, opts)      ! Bisection
    x2 = regfal (lb, ub, fp)            ! Regula Falsi


    print *, "x: ", x1, x2
    print *, "f(x): ", fp(x1), fp(x2)


!   tests for non-enclosing intervals
!   lb = 6.0e-2_real64
!   ub = 7.0e-2_real64
!   x1 = bisect (lb, ub, fp)    passed
!   x2 = regfal (lb, ub, fp)    passed
end program
