!
!   Applied Numerical Analysis                                July 27, 2021
!   ME 2020 FA21
!   Prof. M Diaz-Maldonado
!
!   source: numint.for
!
!   Synopsis:
!   Implements (some) numerical integration techniques.
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

module num_integrators
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none
    private

    interface
        function fun(x) result(f)
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none
            real(kind = real64), intent(in) :: x(:)
            real(kind = real64) :: f(size(x))
        end function
    end interface

    public :: lsum
    public :: rsum
    public :: trap
    contains


        function lsum (lb, ub, n, fp) result(ni)        ! Left Riemann
            procedure(fun), pointer :: fp               ! f(x)
            real(kind = real64), intent(in) :: lb, ub   ! [low, up] limits
            integer(kind = int32), intent(in) :: n      ! intervals
            integer(kind = int32):: i                   ! counter
            real(kind = real64):: ni                    ! numeric integral
            real(kind = real64):: dx                    ! step
            real(kind = real64):: x(n + 1)              ! values

            dx = (ub - lb) / real(n, kind = real64)
            x(:) = [(lb + real(i, kind = real64) * dx, i = 0, n)]
            ni = dx * sum( fp( x(1:n) ) )
            return
        end function


        function rsum (lb, ub, n, fp) result(ni)        ! Right Riemann
            procedure(fun), pointer :: fp
            real(kind = real64), intent(in) :: lb, ub
            integer(kind = int32), intent(in) :: n
            integer(kind = int32):: i
            real(kind = real64):: ni
            real(kind = real64):: dx
            real(kind = real64):: x(n + 1)

            dx = (ub - lb) / real(n, kind = real64)
            x(:) = [(lb + real(i, kind = real64) * dx, i = 0, n)]
            ni = dx * sum( fp( x(2:n+1) ) )
            return
        end function


        function trap (lb, ub, n, fp) result(ni)        ! Trapezoid Method
            procedure(fun), pointer :: fp
            real(kind = real64), intent(in) :: lb, ub
            integer(kind = int32), intent(in) :: n
            integer(kind = int32):: i
            real(kind = real64):: ni
            real(kind = real64):: dx
            real(kind = real64):: x(n + 1)

            dx = (ub - lb) / real(n, kind = real64)
            x(:) = [(lb + real(i, kind = real64) * dx, i = 0, n)]
            ni = sum( fp( x(1:1) ) ) + 2.0_real64 * sum( fp(x(2:n)) ) + &
               & sum( fp( x(n+1:n+1) ) )
            ni = 0.5_real64 * dx * ni
            return
        end function


end module


! Comments:
! Had to resort to some gimmicks to reduce rank-1 arrays into scalars while
! keeping the number of lines of code to a minimum.
!
!
! Warnings:
! Watch out for stack overflows if using too many intervals for numerical
! integration. If you really need that many you would have to allocate the
! arrays in the heap via the allocate statement.
