!
!   author: misael-diaz                                     August 17, 2021
!   source: numint.for
!
!   Synopsis:
!   Implements (some) Ordinary Differential Equation Solvers.
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

module odes
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    private

    interface
        function odefun(t, y) result(f)
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64) :: f
        end function
    end interface

    public :: Euler
    contains


        subroutine Euler (odesol, ti, tf, yi, N, f)
            ! Synopsis: Possible implementation of Euler's explicit method
            integer(kind = int64), intent(in) :: N              ! steps
            real(kind = real64), intent(in) :: ti, tf           ! tspan
            real(kind = real64), intent(in) :: yi               ! y(t = ti)
            real(kind = real64), intent(inout), target :: odesol(N + 1, 2)
            procedure(odefun), pointer :: f                     ! odefun
            integer(kind = int64) :: i                          ! counter
            real(kind = real64) :: dt                           ! time-step
            real(kind = real64), pointer, contiguous :: t(:) => null()
            real(kind = real64), pointer, contiguous :: y(:) => null()

            t => odesol(:, 1)   ! time, t
            y => odesol(:, 2)   ! solution, y(t)
            dt = (tf - ti) / real(N, kind = real64)
            call linspace (t, ti, tf, N + 1)


            y(1) = yi
            do i = 1, N
                y(i + 1) = y(i) + dt * f( t(i), y(i) )
            end do


            return
        end subroutine


        subroutine linspace (x, xi, xf, numel)
            ! Synopsis: Implements a numpy-like linspace method
            integer(kind = int64), intent(in) :: numel
            real(kind = real64), intent(inout) :: x(numel)
            integer(kind = int64) :: i
            real(kind = real64) :: xi, xf, dx

            dx = (xf - xi) / real(numel - 1_int64, kind = real64)

            do i = 1_int64, numel
                x(i) = xi + real(i - 1_int64, kind = real64) * dx
            end do

            return
        end subroutine


end module
