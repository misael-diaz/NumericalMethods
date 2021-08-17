!
!   source: test.for                                        August 17, 2021
!   author: misael-diaz                                 
!
!   Synopsis:
!   Tests the implementation of Ordinary Differential Equation ODE solvers.
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


module odefuns
    ! encapsulates the ODE function f(x) in a module for convenience
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    public
    contains
        function odefun(t, y) result(f)
            ! Synopsis:
            ! Defines the nonlinear function, f(x), to be integrated.
            real(kind = real64), intent(in) :: t
            real(kind = real64), intent(in) :: y
            real(kind = real64), parameter :: k = 1.0_real64
            real(kind = real64):: f

            f = -k * y
            
            return
        end function
end module


program tests
    ! Solves linear Ordinary Differential Equations ODEs numerically.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use odes, only: Euler
    use odefuns, only: odefun
    implicit none

    interface
        subroutine write (odesol, filename)
            use, intrinsic :: iso_fortran_env, only: int64, real64
            implicit none
            integer(kind = int64) :: i, n, unit, iostat
            real(kind = real64), intent(in) :: odesol(:, :)
            character(len=*), intent(in) :: filename
        end subroutine
    end interface

    real(kind = real64):: ti, tf
    real(kind = real64):: yi
    real(kind = real64), allocatable :: odesol(:, :)
    procedure(odefun), pointer :: fp => null()
    integer(kind = int64), parameter :: n = 255
    integer(kind = int64) :: mstat
    character(:), allocatable :: filename


    allocate (odesol(n + 1, 2), stat=mstat)
    if (mstat /= 0) error stop "failed to allocate memory buffers"


    fp => odefun
    ti = 0.0_real64             ! initial time
    tf = 5.0_real64             ! final time
    yi = 1.0_real64             ! initial value, y(t = ti) = yi

    call Euler (odesol, ti, tf, yi, n, fp)

    ! allocates string of characters dynamically
    allocate (character( len = len("output/Euler.dat") ) :: filename, &
        & stat = mstat)
    if (mstat /= 0) error stop "failed to allocate memory for string"

    filename = "output/Euler.dat"
    call write (odesol, filename)

    deallocate (odesol, filename, stat=mstat)
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
