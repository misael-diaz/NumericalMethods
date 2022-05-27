/*
 * Applied Numerical Analysis			         September 20, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: impulse.cpp
 *
 * Synopsis:
 * Obtains the transient response of a first-order Ordinary Differential
 * Equation ODE subject to step and impulse responses:
 *
 * 		y' + k * y = b * u(t),	y(0) = 0,
 * 
 * where k and b are the rate and forcing-constants, respectively, and
 * u(t) is either the unit-step or the unit-impulse function.
 *
 *
 * Copyright (c) 2021 Misael Diaz-Maldonado
 *
 * This file is released under the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 *
 * References:
 * [0] A Gilat and V Subramanian, Numerical Methods for Engineers and
 *     Scientists, 3rd edition.
 * [1] A Koenig and B Moo, Accelerated C++ Practical Programming by
 *     Example.
 * [2] www.cplusplus.com/reference/tuple/tuple
 * [3] www.geeksforgeeks.org/tuples-in-c/
 *
 */

#include <ios>
#include <cmath>
#include <tuple>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <functional>
// defines the rate and forcing constant, respectively
#define RATE 1.0
#define FEXT 1.0
// implements the absolute value function via MACROS
#define absval(x) (signbit(x)? -x: x)
import odes ;

using std::get ;
using std::cout ;
using std::endl ;
using std::tuple ;
using std::string ;
using std::vector ;
using std::signbit ;
using std::ofstream ;
using std::scientific ;
using std::streamsize ;
using std::setprecision ;
using std::runtime_error ;

using ode::iEuler ;
using ode::RK2 ;

// lambda, odefun, RHS of the ODE f(t, y)
std::function < double(double, double) >
odefun = [](double t, double y) -> double {
	const double k = RATE ;
	const double b = FEXT ;
	return (b - k * y) ;
} ;

// analytic expression for the step-response, y(t)
auto fstep = [](double t) -> double {
	const double k  = RATE ;
	const double b  = FEXT ;
	return ( -(b / k) * exp(-k * t) + b / k ) ;
} ;

// analytic expression for the impulse-response, y(t)
auto fimpulse = [](double t) -> double {
	const double k  = RATE ;
	const double b  = FEXT ;
	return  ( b * exp(-k * t) );
} ;

// loggers
void write (const std::tuple <std::vector<double>, std::vector<double>>& );
void write (const std::tuple < std::vector<double>, std::vector<double> >&,
	    const std::string& );

int main() {

	// defines parameters for the solver
	string filename ;
	const int N = 511 ;			// number of time-steps
	const double yi = 0.0 ;			// initial value y(0)
	const double ti = 0.0, tf = 3.0e1 ;	// time span

	// creates placeholders for the numerical solutions
	tuple < vector<double>, vector<double> > odesol_iEuler ;
	tuple < vector<double>, vector<double> > odesol_RK2 ;

	// solves the ODE using the specified method
	odesol_iEuler = iEuler (odesol_iEuler, ti, tf, yi, N, odefun);
	odesol_RK2    = RK2    (odesol_RK2,    ti, tf, yi, N, odefun);
	
	// exports the numerical solutions
	filename = "output/impulse/iEulr.dat" ;
	write (odesol_iEuler, filename);
	filename = "output/impulse/EuRK2.dat" ;
	write (odesol_RK2, filename);
	write (odesol_RK2);
	return 0 ;
}


void write (const tuple < vector<double>, vector<double> >& odesol)
{	// writes numerical results on the output stream

	// unpacks numerical solution from tuple
	const vector<double>& t = get<0>(odesol);
	const vector<double>& y = get<1>(odesol);

	streamsize prec = cout.precision() ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != t.size() ; ++i) {
		cout << scientific << setprecision(15) <<
			t[i] << "\t" << y[i] << "\t" <<
			odefun(t[i], y[i]) << "\t" <<
			absval( (fstep(t[i]) - y[i]) ) << "\t" <<
			absval( (fimpulse(t[i]) - odefun(t[i], y[i])) ) <<
			setprecision(prec) << endl ;
	}
}


void write (const tuple < vector<double>, vector<double> >& odesol,
	    const string& filename)
{	// writes numerical results to the output filestream

	const vector<double>& t = get<0>(odesol);
	const vector<double>& y = get<1>(odesol);

	ofstream out ;
	out.open (filename, std::ios::out);

	if ( !out.is_open() ) {
		throw runtime_error("I/O Error: " + filename);
	}

	streamsize prec = out.precision() ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != t.size() ; ++i) {
		out << scientific << setprecision(15) <<
			t[i] << "\t" << y[i] << "\t" <<
			odefun(t[i], y[i]) << "\t" <<
			absval( (fstep(t[i]) - y[i]) ) << "\t" <<
			absval( (fimpulse(t[i]) - odefun(t[i], y[i])) ) <<
			setprecision(prec) << endl ;
	}

	out.close() ;
}
