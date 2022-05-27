!
!   Applied Numerical Analysis                           September 14, 2021
!   Prof. M Diaz-Maldonado
!
!   Synopsis:
!   Solves a nonlinear Ordinary Differential Equation ODE, which describes
!   the radiation heat transfer from a slab having a uniform temperature
!   distribution:
!
!    cp * rho * L * dT/dt = eps * sigma * (T**4 - T_sur**4),   T(0) = T0,
!
!   where `cp' is the heat capacity, `rho` is the density, `L' is the slab
!   thickness, `T' is the absolute temperature of the slab, `eps' is the
!   emissivity, `sigma' is the Stefan-Boltzmann constant, `T_sur' is the
!   surroundings temperature, `T0' is the initial slab temperature, and
!   `t' the time.
!
!   Note:
!   The absolute temperature (Kelvins) is scaled by a factor of 1e3 and the
!   time is non-dimensionalized by 1.0e9 * eps * sigma / (rho * cp * L) to
!   minimize numerical errors.
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
!   [0] TL Bergman, AS Lavine, FP Incropera, DP DeWitt, Fundamentals
!       of Heat and Mass Transfer, 8th edition.
!   [1] A Gilat and V Subramanian, Numerical Methods for Engineers and
!       Scientists, 3rd edition.
!   [2] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
!   [3] SJ Chapman, FORTRAN for Scientists and Engineers, 4th edition


! MACROS for the surroundings and initial slab temperature, respectively
#define TSUR  0.29315_real64
#define TEMP0 1.2_real64

module radiation_odefuns
    ! defines functions for solving for the transient response
    use, intrinsic :: iso_fortran_env, only: int64, real64
    implicit none
    private
    public :: odefun
    public :: odesol_Write
    contains
        function odefun(t, y, params) result(f)
            ! Synopsis:
            ! RHS of the ODE that describes the slab temperature evolution.
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), intent(in) :: params(:)
            real(kind = real64):: Temp, T_sur
            real(kind = real64):: f

            Temp = y                    ! slab temperature
            T_sur = params(1)           ! surroundings temperature
            f = -(Temp**4 - T_sur**4)   ! RHS, f(t, T)

            return
        end function


        function ferr(t, y) result(f)
            ! Synopsis: Error function based on the analytic solution.
            real(kind = real64):: f, Temp
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), parameter :: T0    = TEMP0
            real(kind = real64), parameter :: T_sur = TSUR
            real(kind = real64), parameter :: C     = 1.0_real64 / &
                & (4.0_real64 * T_sur**3)

            Temp = y
            f = ( C * (2.0 * ( datan(Temp / T_sur) - datan(T0 / T_sur) ) +&
                & dlog( (T_sur + Temp) / (T_sur + T0) ) - &
                & dlog( (T_sur - Temp) / (T_sur - T0) ) ) - t )

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
                err = dabs( ferr( t(i), y(i) ) )
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
    use radiation_odefuns, only: odefun
    use radiation_odefuns, only: write => odesol_Write
    implicit none


    real(kind = real64):: T_sur
    real(kind = real64):: yi                                    ! y(t = ti)
    real(kind = real64):: ti, tf                                ! tspan
    integer(kind = int64), parameter :: n = 65535               ! num steps
    real(kind = real64), allocatable, target :: odesol(:, :, :) ! solution
    type(ODE_solverParams), allocatable, target :: params       ! ODE Param
    real(kind = real64), pointer, contiguous :: p_odesol(:, :) => null()
    procedure(odefun), pointer :: fp => null()
    integer(kind = int64) :: mstat
    character(:), allocatable :: filename
    character(*), parameter :: BaseFilename = &
        & "output/heat-transfer/radiation/Euler.dat"


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
    tf = 5.0e1_real64           ! final time
    yi = TEMP0                  ! initial value, y(t = ti) = yi
    T_sur = TSUR

    associate (prms => params % prms)
        prms(1) = T_sur
    end associate


    !! solves the ODEs with the specified method and exports results
    p_odesol => odesol(:, :, 1)
    call iEuler ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/heat-transfer/radiation/iEulr.dat"
    call write (p_odesol, filename)

    p_odesol => odesol(:, :, 2)
    call RK2    ( p_odesol, ti, tf, yi, n, fp, c_loc(params) )
    filename = "output/heat-transfer/radiation/EuRK2.dat"
    call write (p_odesol, filename)


    !! frees memory buffers
    deallocate (params, odesol, filename, stat=mstat)
    if (mstat /= 0) error stop "unexpected memory deallocation error"
end program
