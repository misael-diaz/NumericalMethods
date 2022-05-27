/*
 * Applied Numerical Analysis			           October 02, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: isothermal-second-order-kinetics-batch-reactor.cpp
 *
 * Synopsis:
 * Obtains the transient response of a nonlinear Ordinary Differential
 * Equation ODE, which describes the depletion of a chemical species by a
 * chemical reaction of second-order kinetics carried out in a isothermal
 * batch reactor:
 *
 *                      y' = -beta * y**2,	y(0) = 1,
 *
 * where `beta' is the effective reaction rate constant, beta = k * Ca0,
 * where `k' is the reaction rate constant, `Ca0' is the initial reactant
 * concentration (moles / volume), `y' is the non-dimensional reactant
 * concentration y(t) = Ca(t) / Ca0, and `t' is the time.
 *
 * It's assumed that the chemical reaction takes place in a liquid so that
 * the volume of the reactor can be regarded as constant.
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
 * [2] HS Fogler, Elements of Chemical Reaction Engineering, 6th edition
 * [3] CA Kluever, Dynamic Systems: Modeling, Simulation, and Control
 * [4] www.cplusplus.com/reference/tuple/tuple
 * [5] www.geeksforgeeks.org/tuples-in-c/
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
#define BETA 1.0
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
	const double beta = BETA ;
	return (-beta * y * y);
} ;

// analytic expression for the non-dimensional reactant concentration, y(t)
auto fsol = [](double t) -> double {
	return ( 1.0 / (1.0 + BETA * t) );
} ;

// loggers
void write (const std::tuple <std::vector<double>, std::vector<double>>& );
void write (const std::tuple < std::vector<double>, std::vector<double> >&,
	    const std::string& );

int main() {

	// defines parameters for the solver
	string filename ;
	const int N = 255 ;			// number of time-steps
	const double yi = 1.0 ;			// initial value y(0)
	const double ti = 0.0, tf = 10.0 ;	// time span

	// creates placeholders for the numerical solutions
	tuple < vector<double>, vector<double> > odesol_iEuler ;
	tuple < vector<double>, vector<double> > odesol_RK2 ;

	// solves the ODE using the specified method
	odesol_iEuler = iEuler (odesol_iEuler, ti, tf, yi, N, odefun);
	odesol_RK2    = RK2    (odesol_RK2,    ti, tf, yi, N, odefun);

	// exports the numerical solutions
	filename = "output/cheme/kinetics/iEulr.dat" ;
	write (odesol_iEuler, filename);
	filename = "output/cheme/kinetics/EuRK2.dat" ;
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
			absval( (fsol(t[i]) - y[i]) ) << "\t" <<
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
			absval( (fsol(t[i]) - y[i]) ) << "\t" <<
			setprecision(prec) << endl ;
	}

	out.close() ;
}
