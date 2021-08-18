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
!
!
!   References:
!   [0] SJ Chapman, FORTRAN for Scientists and Engineers, 4th edition.
!   [1] A Gilat and V Subramanian, Numerical Methods for Engineers and
!       Scientists, 3rd edition.
!   [2] gcc.gnu.org/onlinedocs/gfortran/Working-with-Pointers.html

module odes
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use nlsolvers, only: fzero
    implicit none
    private

    interface
        function odefun(t, y, params) result(f)
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), intent(in) :: params(:)
            real(kind = real64) :: f
        end function

        function objfun (yn) result(objf)
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none
            real(kind = real64), intent(in) :: yn
            real(kind = real64):: objf
        end function
    end interface

    type, public :: ODE_solverParams
        real(kind = real64), allocatable :: prms(:)
        contains
            final :: finalizer
    end type    ! ODE Solver Parameters

    public :: Euler
    public :: iEuler
    contains


        subroutine Euler (odesol, ti, tf, yi, N, f, vprms)
            ! Synopsis: Possible implementation of Euler's explicit method
            real(kind = real64), intent(in) :: ti, tf           ! tspan
            real(kind = real64), intent(in) :: yi               ! y(t = ti)
            integer(kind = int64), intent(in) :: N              ! steps
            procedure(odefun), pointer :: f                     ! odefun
            type(c_ptr), intent(in), value :: vprms             ! void*
            integer(kind = int64) :: i                          ! counter
            real(kind = real64) :: dt                           ! time-step
            real(kind = real64), intent(inout), target :: odesol(N + 1, 2)
            type(ODE_solverParams), pointer :: params => null()
            real(kind = real64), pointer, contiguous :: t(:) => null()
            real(kind = real64), pointer, contiguous :: y(:) => null()
            real(kind = real64), pointer, contiguous :: prms(:) => null()

            call c_f_pointer (vprms, params)    ! binds ode params to void*
            prms => params % prms               ! binds to the param array

            t => odesol(:, 1)   ! time, t
            y => odesol(:, 2)   ! solution, y(t)
            dt = (tf - ti) / real(N, kind = real64)
            call linspace (t, ti, tf, N + 1)


            y(1) = yi
            do i = 1, N
                y(i + 1) = y(i) + dt * f( t(i), y(i), prms )
            end do

            return
        end subroutine


        subroutine iEuler (odesol, ti, tf, yi, N, fp_odefun, vprms)
            ! Synopsis: Possible implementation of Euler's implicit method
            real(kind = real64), intent(in) :: ti, tf           ! tspan
            real(kind = real64), intent(in) :: yi               ! y(t = ti)
            integer(kind = int64), intent(in) :: N              ! steps
            procedure(odefun), intent(in), pointer :: fp_odefun ! odefun
            procedure(objfun), pointer :: fp_objfun => null()   ! objfun
            type(c_ptr), intent(in), value :: vprms             ! void*
            integer(kind = int64) :: i                          ! counter
            real(kind = real64) :: dt                           ! time-step
            real(kind = real64) :: K1, K2                       ! slopes
            real(kind = real64) :: y_lb, y_ub                   ! bounds
            real(kind = real64), intent(inout), target :: odesol(N + 1, 2)
            type(ODE_solverParams), pointer :: params => null()
            real(kind = real64), pointer, contiguous :: t(:) => null()
            real(kind = real64), pointer, contiguous :: y(:) => null()
            real(kind = real64), pointer, contiguous :: prms(:) => null()

            fp_objfun => objfun
            call c_f_pointer (vprms, params)
            prms => params % prms

            t => odesol(:, 1)
            y => odesol(:, 2)
            dt = (tf - ti) / real(N, kind = real64)
            call linspace (t, ti, tf, N + 1)


            y(1) = yi
            do i = 1, N
                ! bounds the solution y(i+1) from below and above
                K1 = fp_odefun ( t(i), y(i), prms )
                K2 = fp_odefun ( t(i) + dt, y(i) + K1 * dt, prms )
                y_lb = y(i) + dt * K1
                y_ub = y(i) + dt * K2
                ! solves for y(i+1) iteratively
                y(i + 1) = fzero (y_lb, y_ub, fp_objfun)
            end do


            return
            contains
                function objfun (yn) result(objf)
                    ! objective function supplied to nonlinear solver
                    real(kind = real64), intent(in) :: yn
                    real(kind = real64):: objf

                    objf = yn - y(i) - dt * fp_odefun( t(i+1), yn, prms )

                    return
                end function
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



        subroutine finalizer (params)
            ! Destroys the fields of objects of type ODE Solver Parameters.
            type(ODE_solverParams), intent(inout) :: params
            integer(kind = int64):: mstat

            if ( allocated(params % prms) ) then
                deallocate(params % prms, stat=mstat)
                if (mstat /= 0) error stop "finalizer: unexpected error"
            end if

            return
        end subroutine


end module


! Comments:
! The type(c_ptr) is used as a FORTRAN analogue of C's universal pointer
! (void *). Note that it's passed by value to the numerical techniques
! so that it behaves just like void* (otherwise it behaves as void**,
! see GCC reference above). To achieve that it's also necessary to use
! c_loc on the caller side as shown in the main program.
!
! The rationale for using it is that it allows users to apply the numerical
! techniques with their own derived-types. The odefun subroutine just needs
! to "unpack" the needed parameters to evaluate the expression f(t, y).
!
! Defining a finalizer saves the user from manually deallocating the fields
! of derived-types.
