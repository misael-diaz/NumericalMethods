!
!   Applied Numerical Analysis                              August 28, 2021
!   Prof. M Diaz-Maldonado
!
!   Synopsis:
!   Solves for the transient response of a first-order dynamic system
!   subject to a unit-ramp input:
!
!                       y' + k * y = b * u(t),
!
!   where k is the rate, b is the forcing constant, and u(t) = t is the
!   unit-ramp function for t > 0.
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
!   [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
!       Scientists, 3rd edition.
!   [1] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
!   [2] SJ Chapman, FORTRAN for Scientists and Engineers, 4th edition

! MACROS for the rate and forcing constants, and the initial value y(0)
#define RATE 1.0_real64
#define FEXT 1.0_real64
#define YINI 0.0_real64

module ramp_odefuns
    ! defines functions for solving for the ramp response
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: odefun
    public :: odesol_Write
    contains
        function odefun(t, y, params) result(f)
            ! Synopsis:
            ! RHS of the first-order ODE subject to a unit-ramp input.
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), intent(in) :: params(:)
            real(kind = real64):: k, b
            real(kind = real64):: f

            k = params(1)
            b = params(2)
            f = b * t - k * y
            
            return
        end function


        function framp(t) result(y)
            ! Synopsis: Analytic expression for the ramp response y(t).
            real(kind = real64) :: y
            real(kind = real64), intent(in) :: t
            real(kind = real64), parameter :: k  = RATE
            real(kind = real64), parameter :: b  = FEXT

            y = (b / (k * k) * (dexp(-k * t) - 1.0_real64) + b / k * t)
            
            return
        end function


        subroutine odesol_Write (odesol, filename)
            ! Synopsis: Writes the numerical solution (t, y) to a data file
            integer(kind = int64) :: i, nrows, unit, iostat
            real(kind = real64) :: err
            real(kind = real64), intent(in), target :: odesol(:, :)
            real(kind = real64), pointer, contiguous :: t(:) => null()
            real(kind = real64), pointer, contiguous :: y(:) => null()
            character(len=*), intent(in) :: filename

            unit = 100
            open(unit=unit, file=filename, action='write', iostat = iostat)
            if (iostat /= 0) error stop ("I/O error: " // filename)

            t => odesol(:, 1)
            y => odesol(:, 2)

            nrows = size (array = odesol, dim = 1, kind = int64)
            do i = 1, nrows
                err = dabs( framp( t(i) ) - y(i) )
                write (unit, '(3E25.16)') t(i), y(i), err
            end do

            close(unit)
            return
        end subroutine
end module


program tests
    ! Solves first-order Ordinary Differential Equations ODEs numerically.
    use, intrinsic :: iso_c_binding, only: c_loc
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use odes, only: iEuler, RK2, ODE_solverParams
    use ramp_odefuns, only: odefun
    use ramp_odefuns, only: write => odesol_Write
    implicit none


    real(kind = real64):: k                                     ! rate
    real(kind = real64):: b                                     ! forcing
    real(kind = real64):: yi                                    ! y(t = ti)
    real(kind = real64):: ti, tf                                ! tspan
    integer(kind = int64), parameter :: n = 255                 ! num steps
    real(kind = real64), allocatable, target :: odesol(:, :, :) ! solution
    type(ODE_solverParams), allocatable, target :: params       ! ODE Param
    real(kind = real64), pointer, contiguous :: p_odesol(:, :) => null()
    procedure(odefun), pointer :: fp => null()
    integer(kind = int64) :: mstat
    character(:), allocatable :: filename
    character(*), parameter :: BaseFilename = "output/ramp/Euler.dat"


    !! memory allocations
    allocate (params, odesol(n + 1, 2, 2), stat=mstat)
    if (mstat /= 0) error stop "failed to allocate memory buffers"

    allocate (params % prms(2), stat=mstat)
    if (mstat /= 0) error stop "failed to allocate array of parameters"

    allocate (character( len = len(BaseFilename) ) :: filename, &
        & stat = mstat) ! dynamically allocated string of characters
    if (mstat /= 0) error stop "failed to allocate memory for string"


    !! initializations
    fp => odefun
    ti = 0.0_real64             ! initial time
    tf = 5.0_real64             ! final time
    yi = YINI                   ! initial value, y(t = ti) = yi
    k  = RATE                   ! rate constant
    b  = FEXT                   ! forcing constant

    associate (prms => params % prms)
        prms(1) = k
        prms(2) = b
    end associate


    !! solves the ODEs with the specified method and exports results
    p_odesol => odesol(:, :, 1)
    call iEuler ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/ramp/iEulr.dat"
    call write (p_odesol, filename)

    p_odesol => odesol(:, :, 2)
    call RK2    ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/ramp/EuRK2.dat"
    call write (p_odesol, filename)


    !! frees memory buffers
    deallocate (params, odesol, filename, stat=mstat)
    if (mstat /= 0) error stop "unexpected memory deallocation error"
end program


subroutine write (odesol, filename)
    ! Synopsis: Writes the numerical solution (t, y) to a data file
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    integer(kind = int64) :: i, nrows, unit, iostat
    real(kind = real64), intent(in) :: odesol(:, :)
    character(len=*), intent(in) :: filename

    unit = 100
    open (unit=unit, file=filename, action='write', iostat = iostat)
    if (iostat /= 0) error stop ("I/O error: " // filename)

    nrows = size (array = odesol, dim = 1, kind = int64)
    do i = 1, nrows
        write (unit, '(2E25.15)') odesol(i, :)
    end do

    close(unit)
    return
end subroutine
