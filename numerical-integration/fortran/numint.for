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
            real(kind = real64), intent(in) :: x
            real(kind = real64) :: f
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
            real(kind = real64):: x                     ! values
            real(kind = real64):: s                     ! accumulator

            dx = (ub - lb) / real(n, kind = real64)

            s = 0.0_real64
            do i = 0, n - 1
                x = lb + real(i, kind = real64) * dx
                s = s + fp(x)
            end do
            ni = dx * s

            return
        end function


        function rsum (lb, ub, n, fp) result(ni)        ! Right Riemann
            procedure(fun), pointer :: fp
            real(kind = real64), intent(in) :: lb, ub
            integer(kind = int32), intent(in) :: n
            integer(kind = int32):: i
            real(kind = real64):: ni
            real(kind = real64):: dx
            real(kind = real64):: x
            real(kind = real64):: s

            dx = (ub - lb) / real(n, kind = real64)

            s = 0.0_real64
            do i = 1, n
                x = lb + real(i, kind = real64) * dx
                s = s + fp(x)
            end do
            ni = dx * s

            return
        end function


        function trap (lb, ub, n, fp) result(ni)        ! Trapezoid Method
            procedure(fun), pointer :: fp
            real(kind = real64), intent(in) :: lb, ub
            integer(kind = int32), intent(in) :: n
            integer(kind = int32):: i
            real(kind = real64):: ni
            real(kind = real64):: dx
            real(kind = real64):: x
            real(kind = real64):: s

            dx = (ub - lb) / real(n, kind = real64)

            s = 0.0_real64
            do i = 1, n - 1
                x = lb + real(i, kind = real64) * dx
                s = s + 2.0_real64 * fp(x)
            end do
            ni = 0.5_real64 * dx * ( fp(lb) + s + fp(ub) )

            return
        end function


end module
