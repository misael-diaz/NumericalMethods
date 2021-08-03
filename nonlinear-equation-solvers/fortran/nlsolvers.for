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

    type, public :: nls_conf    !! conf[iguration] struct
        real(kind = real64) :: tol = TOL
        integer(kind = int32) :: max_iter = MAX_ITER
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
    contains


        function bisect (lb, ub, fp, opts) result(x)    ! Bisection
            real(kind = real64), intent(in) :: lb, ub   ! [low, up] bounds
            procedure(fun), pointer :: fp               ! f(x)
            real(kind = real64):: a, b                  ! bounds aliases
            real(kind = real64):: x                     ! root approximate
            integer(kind = int32):: n, maxit            ! count && max iter
            real(kind = real64):: t                     ! tolerance
            type(nls_conf), intent(in), optional :: opts

            call bracket_check (lb, ub, fp)

            if ( present(opts) ) then
                t     = opts % tol
                maxit = opts % max_iter
            else
                t     = TOL
                maxit = MAX_ITER
            end if

            ! checks bounds
            if (lb < ub) then
                a = lb
                b = ub
            else
                a = ub
                b = lb
            end if

            n = 1
            x = 0.5_real64 * (a + b)
            do while ( n /= maxit .and. abs( fp(x) ) > t )

                ! selects the bracketing interval
                if ( fp(a) * fp(x) < 0.0_real64 ) then
                    b = x
                else
                    a = x
                end if

                x = 0.5_real64 * (a + b)
                n = n + 1
            end do

            call report (n)

            return
        end function


        function regfal (lb, ub, fp) result(x)          ! Regula Falsi
            real(kind = real64), intent(in) :: lb, ub   ! [low, up] bounds
            procedure(fun), pointer :: fp               ! f(x)
            real(kind = real64):: a, b                  ! bounds aliases
            real(kind = real64):: x                     ! root approximate
            integer(kind = int32):: n                   ! iterations

            call bracket_check (lb, ub, fp)

            ! checks bounds
            if (lb < ub) then
                a = lb
                b = ub
            else
                a = ub
                b = lb
            end if

            n = 1
            x = ( a * fp(b) - b * fp(a) ) / ( fp(b) - fp(a) )
            do while ( n /= MAX_ITER .and. abs( fp(x) ) > TOL )

                ! selects the bracketing interval
                if ( fp(a) * fp(x) < 0.0_real64 ) then
                    b = x
                else
                    a = x
                end if

                x = ( a * fp(b) - b * fp(a) ) / ( fp(b) - fp(a) )
                n = n + 1
            end do

            call report (n)

            return
        end function


        subroutine report(n)
            integer(kind = int32), intent(in) :: n

            if (n /= MAX_ITER) then
                print *, "solution found in ", n, " iterations"
            else
                print *, "maximum number of iterations has been " // &
                       & "reached, try again with a narrower interval"
            end if

            return
        end subroutine


        subroutine bracket_check (lb, ub, fp)
            ! complains if there's no root in given interval [lb, ub].
            real(kind = real64), intent(in) :: lb, ub
            procedure(fun), pointer :: fp
            character(len=*), parameter :: errmsg = &
                & "no root exists in given interval"

            if ( fp(lb) * fp(ub) > 0.0_real64 ) then
                error stop errmsg
            end if

            return
        end subroutine


end module



! TODO:
! [x] guard against non-existing root in given interval using error stop
! [ ] consider adding tutorial-like comments at the end of the source
