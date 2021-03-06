/*
 * Applied Numerical Analysis			           October 02, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: sinusoid.cpp
 *
 * Synopsis:
 * Obtains the transient response of a first-order Ordinary Differential
 * Equation ODE subject to a unit-sinusoid input:
 *
 * 		y' + k * y = b * u(t),	y(0) = 0,
 *
 * where k and b are the rate and forcing-constants, respectively, and
 * u(t) = sin(w * t) is the unit-sinusoid input function of frequency `w`.
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
// defines the rate and forcing constant, and frequency, respectively
#define RATE  1.0
#define FEXT  1.0
#define OMEGA 1.0
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
	const double w = OMEGA ;
	return (b * sin(w * t) - k * y) ;
} ;

// analytic expression for the sinusoid-response, y(t)
auto fsinusoid = [](double t) -> double {
	const double k  = RATE ;
	const double b  = FEXT ;
	const double w  = OMEGA ;
        // amplitudes of the sinusoidal response
        double A0 = -w * b / (w * w + k * k);
        double A1 =  k * b / (w * w + k * k);
        return ( A0 * ( cos(w * t) - exp(-k * t) ) + A1 * sin(w * t) );
} ;

// loggers
void write (const std::tuple <std::vector<double>, std::vector<double>>& );
void write (const std::tuple < std::vector<double>, std::vector<double> >&,
	    const std::string& );

int main() {

	// defines parameters for the solver
	string filename ;
	const int N = 255 ;			// number of time-steps
	const double yi = 0.0 ;			// initial value y(0)
	const double ti = 0.0, tf = 5.0 ;	// time span

	// creates placeholders for the numerical solutions
	tuple < vector<double>, vector<double> > odesol_iEuler ;
	tuple < vector<double>, vector<double> > odesol_RK2 ;

	// solves the ODE using the specified method
	odesol_iEuler = iEuler (odesol_iEuler, ti, tf, yi, N, odefun);
	odesol_RK2    = RK2    (odesol_RK2,    ti, tf, yi, N, odefun);

	// exports the numerical solutions
	filename = "output/sinusoid/iEulr.dat" ;
	write (odesol_iEuler, filename);
	filename = "output/sinusoid/EuRK2.dat" ;
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
			absval( (fsinusoid(t[i]) - y[i]) ) << "\t" <<
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
			absval( (fsinusoid(t[i]) - y[i]) ) << "\t" <<
			setprecision(prec) << endl ;
	}

	out.close() ;
}
