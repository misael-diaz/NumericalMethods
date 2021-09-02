!
!   Applied Numerical Analysis                            September 02, 2021
!   Prof. M Diaz-Maldonado
!
!
!   Synopsis:
!   Obtains the transient response of a nonlinear Ordinary Differential
!   Equation ODE, which describes the depletion of a chemical species by a
!   chemical reaction of second-order kinetics carried out in a isothermal
!   batch reactor:
!
!                             y' = -beta * y**2,
!
!   where `beta' is the effective reaction rate constant, beta = k * Ca0,
!   where `k' is the reaction rate constant, `Ca0' is the initial reactant
!   concentration (moles / volume), `y' is the non-dimensional reactant
!   concentration y(t) = Ca(t) / Ca0, and `t' is the time.
!
!   It's assumed that the chemical reaction takes place in a liquid so that
!   the volume of the reactor can be regarded as constant.
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
!   [2] HS Fogler, Elements of Chemical Reaction Engineering, 6th edition
!   [3] SJ Chapman, FORTRAN for Scientists and Engineers, 4th edition


! MACROS for the effective reaction rate and the initial value y(0) = 1
#define BETA  1.0_real64
#define YINI  1.0_real64


module kinetics_odefuns
    ! defines functions to obtain the reactant concentration y(t)
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: odefun
    public :: odesol_Write
    contains
        function odefun(t, y, params) result(f)
            ! Synopsis: RHS of the nonlinear ODE (molar balance equation).
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), intent(in) :: params(:)
            real(kind = real64):: beta
            real(kind = real64):: f

            beta = params(1)
            f = -beta * y**2
            
            return
        end function


        function fsol(t) result(y)
            ! Synopsis: The analytic expression for the concentration y(t).
            real(kind = real64) :: y
            real(kind = real64), intent(in) :: t
            real(kind = real64), parameter :: beta  = BETA
                
            y = (  1.0_real64 / (1.0_real64 + beta * t)  )
            
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
                err = dabs( fsol( t(i) ) - y(i) )
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
    use kinetics_odefuns, only: odefun
    use kinetics_odefuns, only: write => odesol_Write
    implicit none


    real(kind = real64):: beta                                  ! rate
    real(kind = real64):: yi                                    ! y(t = ti)
    real(kind = real64):: ti, tf                                ! tspan
    integer(kind = int64), parameter :: n = 255                 ! num steps
    real(kind = real64), allocatable, target :: odesol(:, :, :) ! solution
    type(ODE_solverParams), allocatable, target :: params       ! ODE Param
    real(kind = real64), pointer, contiguous :: p_odesol(:, :) => null()
    procedure(odefun), pointer :: fp => null()
    integer(kind = int64) :: mstat
    character(:), allocatable :: filename
    character(*), parameter :: BaseFilename = &
        & "output/cheme/kinetics/Euler.dat"


    !! memory allocations
    allocate (params, odesol(n + 1, 2, 2), stat=mstat)
    if (mstat /= 0) error stop "failed to allocate memory buffers"

    allocate (params % prms(1), stat=mstat)
    if (mstat /= 0) error stop "failed to allocate array of parameters"

    allocate (character( len = len(BaseFilename) ) :: filename, &
        & stat = mstat) ! dynamically allocated string of characters
    if (mstat /= 0) error stop "failed to allocate memory for string"


    !! initializations
    fp => odefun
    ti = 0.0_real64             ! initial time
    tf = 10.0_real64            ! final time
    yi = YINI                   ! initial value, y(t = ti) = yi
    beta = BETA                 ! effective reaction rate constant

    associate (prms => params % prms)
        prms(1) = beta
    end associate


    !! solves the ODEs with the specified method and exports results
    p_odesol => odesol(:, :, 1)
    call iEuler ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/cheme/kinetics/iEulr.dat"
    call write (p_odesol, filename)

    p_odesol => odesol(:, :, 2)
    call RK2    ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/cheme/kinetics/EuRK2.dat"
    call write (p_odesol, filename)


    !! frees memory buffers
    deallocate (params, odesol, filename, stat=mstat)
    if (mstat /= 0) error stop "unexpected memory deallocation error"

end program
