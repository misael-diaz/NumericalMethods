!
!   Applied Numerical Analysis                                July 26, 2021
!   ME 2020 FA21
!   Prof. M Diaz-Maldonado
!
!   source: nlsolvers.for
!
!
!   Synopsis:
!   Implements (some) nonlinear solvers.
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

module nlsolvers
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    ! constants (defaults for the tolerance and max number of iterations)
    real(kind = real64), parameter :: TOL = 1.0e-8_real64
    integer(kind = int32), parameter :: MAX_ITER = 100
    logical(kind = int32), parameter :: VERBOSE  = .false.

    type, public :: nls_conf    !! conf[iguration] struct
        real(kind = real64) :: tol = TOL
        integer(kind = int32) :: max_iter = MAX_ITER
        logical(kind = int32) :: verbose  = VERBOSE
    end type

    interface
        function fun(x) result(f)
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none
            real(kind = real64), intent(in) :: x
            real(kind = real64) :: f
        end function
    end interface

    public :: bisect
    public :: regfal
    public :: shifter
    contains


        function bisect (lb, ub, fp, opts) result(x)    ! Bisection
            real(kind = real64), intent(in) :: lb, ub   ! [low, up] bounds
            procedure(fun), pointer :: fp               ! f(x)
            real(kind = real64):: a, b                  ! bounds aliases
            real(kind = real64):: x                     ! root approximate
            integer(kind = int32):: n                   ! count
            type(nls_conf) :: conf
            type(nls_conf), intent(in), optional :: opts
            character(len=*), parameter :: nm = "Bisection"

            call bracket_check (lb, ub, fp, nm)
            call bounds_check  (lb, ub, a, b)
            call optset(conf, opts)

            n = 1
            x = 0.5_real64 * (a + b)
            do while (n /= conf % max_iter .and. abs( fp(x) ) > conf % tol)
                call bisector (a, b, x, fp)
                n = n + 1
            end do

            call report (n, nm, conf)

            return
        end function


        subroutine bisector (lb, ub, x, fp)
            ! Synopsis:
            ! Bisects the bracketing interval [lb, ub] and returns
            ! its middle value as an estimate of the root of the
            ! nonlinear function f(x). As a side-effect the function
            ! updates the bounds of the bracketing interval.
            real(kind = real64), intent(inout) :: lb, ub
            real(kind = real64), intent(inout) :: x
            procedure(fun), pointer :: fp

            if ( fp(lb) * fp(x) < 0.0_real64 ) then
                ub = x
            else
                lb = x
            end if

            x = 0.5_real64 * (lb + ub)
            return
        end subroutine


        function regfal (lb, ub, fp, opts) result(x)    ! Regula Falsi
            real(kind = real64), intent(in) :: lb, ub
            procedure(fun), pointer :: fp
            real(kind = real64):: a, b
            real(kind = real64):: x
            integer(kind = int32):: n
            type(nls_conf) :: conf
            type(nls_conf), intent(in), optional :: opts
            character(len=*), parameter :: nm = "Regula Falsi"

            call bracket_check (lb, ub, fp, nm)
            call bounds_check  (lb, ub, a, b)
            call optset(conf, opts)

            n = 1
            x = ( a * fp(b) - b * fp(a) ) / ( fp(b) - fp(a) )
            do while (n /= conf % max_iter .and. abs( fp(x) ) > conf % tol)
                call interp (a, b, x, fp)
                n = n + 1
            end do

            call report (n, nm, conf)

            return
        end function


        subroutine interp (lb, ub, x, fp)
            ! Synopsis:
            ! As bisector but estimates the root via linear interpolation.
            real(kind = real64), intent(inout) :: lb, ub
            real(kind = real64), intent(inout) :: x
            procedure(fun), pointer :: fp

            if ( fp(lb) * fp(x) < 0.0_real64 ) then
                ub = x
            else
                lb = x
            end if

            x = ( lb * fp(ub) - ub * fp(lb) ) / ( fp(ub) - fp(lb) )
            return
        end subroutine


        function shifter (lb, ub, fp, opts) result(x)   ! Shifter Method
            real(kind = real64), intent(in) :: lb, ub
            procedure(fun), pointer :: fp
            real(kind = real64):: a, b
            real(kind = real64):: x1, x2, x
            integer(kind = int32):: n
            type(nls_conf) :: conf
            type(nls_conf), intent(in), optional :: opts
            character(len=*), parameter :: nm = "Shifter"

            call bracket_check (lb, ub, fp, nm)
            call bounds_check  (lb, ub, a, b)
            call optset(conf, opts)

            n = 1
            x1 = 0.5_real64 * (a + b)
            x2 = ( a * fp(b) - b * fp(a) ) / ( fp(b) - fp(a) )
            ! selects the approximate (presumably) closer to the root
            if ( abs(fp(x1)) < abs(fp(x2)) ) then
                x = x1
            else
                x = x2
            end if

            do while (n /= conf % max_iter .and. abs( fp(x) ) > conf % tol)
                call shift (a, b, x, fp)
                n = n + 1
            end do

            call report (n, nm, conf)

            return
        end function


        subroutine shift (lb, ub, x, fp)
            ! Synopsis:
            ! As bisector but shifts towards bisection or interpolation
            ! depending on which yields an approximate closer to the root.
            real(kind = real64), intent(inout) :: lb, ub
            real(kind = real64), intent(inout) :: x
            real(kind = real64) :: x1, x2
            procedure(fun), pointer :: fp

            if ( fp(lb) * fp(x) < 0.0_real64 ) then
                ub = x
            else
                lb = x
            end if

            x1 = 0.5_real64 * (lb + ub)
            x2 = ( lb * fp(ub) - ub * fp(lb) ) / ( fp(ub) - fp(lb) )

            if ( abs(fp(x1)) < abs(fp(x2)) ) then
                x = x1
            else
                x = x2
            end if

            return
        end subroutine


        subroutine report (n, name, conf)
            type(nls_conf), intent(in) :: conf
            integer(kind = int32), intent(in) :: n
            character(len=*), intent(in) :: name
            character(len=*), parameter :: errmsg = &
                & "method needs more iterations for convergence"

            if (n /= conf % max_iter) then
                if (conf % verbose) then
                    print *, name // " Method: "
                    print *, "solution found in ", n, " iterations"
                end if
            else
                error stop (name // " " // errmsg)
            end if

            return
        end subroutine


        subroutine bracket_check (lb, ub, fp, name)
            ! complains if there's no root in given interval [lb, ub].
            real(kind = real64), intent(in) :: lb, ub
            procedure(fun), pointer :: fp
            character(len=*), intent(in) :: name
            character(len=*), parameter :: errmsg = &
                & "No root exists in given interval"

            if ( fp(lb) * fp(ub) > 0.0_real64 ) then
                error stop (name // " " // new_line('N') // errmsg)
            end if

            return
        end subroutine


        subroutine optset (conf, opts)
            ! Synopsis:
            ! Sets the solver conf[iguration] options.
            type(nls_conf), intent(in), optional :: opts
            type(nls_conf), intent(out) :: conf

            if ( present(opts) ) then
                conf % tol      = opts % tol
                conf % max_iter = opts % max_iter
                conf % verbose  = opts % verbose
            else
                conf % tol      = TOL
                conf % max_iter = MAX_ITER
                conf % verbose  = VERBOSE
            end if

            return
        end subroutine


        subroutine bounds_check (lb, ub, a, b)
            ! Synopsis: Ensures the lower bound is less than the upper one.
            real(kind = real64), intent(in) :: lb, ub
            real(kind = real64), intent(out) :: a, b

            if (lb < ub) then
                a = lb
                b = ub
            else
                a = ub
                b = lb
            end if

            return
        end subroutine


end module



! TODO:
! [x] guard against non-existing root in given interval using error stop
! [ ] consider adding tutorial-like comments at the end of the source
