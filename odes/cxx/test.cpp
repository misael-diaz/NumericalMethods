/*
 * Applied Numerical Analysis			            August 10, 2021
 * ME 2020 FA21
 * Prof. M Diaz-Maldonado
 *
 * source: test.cpp
 *
 * Synopsis:
 * Tests (some) Ordinary Differential Equation ODE Solvers.
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
import odes ;

using std::get ;
using std::cout ;
using std::endl ;
using std::tuple ;
using std::string ;
using std::vector ;
using std::ofstream ;
using std::scientific ;
using std::streamsize ;
using std::setprecision ;
using std::runtime_error ;

using ode::Euler ;
using ode::iEuler ;
using ode::RK2 ;


void write (const std::tuple <std::vector<double>, std::vector<double>>& );
void write (const std::tuple < std::vector<double>, std::vector<double> >&,
	    const std::string& );

int main() {

	const int N = 255 ;
	const double yi = 1.0 ;
	const double ti = 0.0, tf = 5.0 ;
	tuple < vector<double>, vector<double> > odesol_Euler ;
	tuple < vector<double>, vector<double> > odesol_iEuler ;
	tuple < vector<double>, vector<double> > odesol_RK2 ;

	// lambda, odefun, RHS of the ODE f(t, y)
	std::function < double(double, double) >
	f = [](double t, double y) -> double {
		const double k = 1.0 ;
		return (-k * y) ;
       	} ;

	// solves the ODE using the specified method
	odesol_Euler  =  Euler(odesol_Euler,  ti, tf, yi, N, f) ;
	odesol_iEuler = iEuler(odesol_iEuler, ti, tf, yi, N, f) ;
	odesol_RK2    = RK2   (odesol_RK2,    ti, tf, yi, N, f) ;
	
	// exports numerical solutions
	string filename = "output/Euler.dat" ;
	write (odesol_Euler, filename) ;

	filename = "output/iEulr.dat" ;
	write (odesol_iEuler, filename) ;

	filename = "output/EuRK2.dat" ;
	write (odesol_RK2, filename) ;

	write (odesol_iEuler) ;
	return 0 ;
}


void write (const tuple < vector<double>, vector<double> >& odesol)
{	// writes numerical results on the output stream

	// unpacks numerical solution in tuple
	const vector<double>& t = get<0>(odesol) ;
	const vector<double>& y = get<1>(odesol) ;

	streamsize prec = cout.precision() ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != t.size() ; ++i) {
		cout << scientific << setprecision(15) << t[i] << "\t" 
		     << y[i] << setprecision(prec) << endl ;
	}
}


void write (const tuple < vector<double>, vector<double> >& odesol,
	    const string& filename)
{	// writes numerical results to the output filestream

	const vector<double>& t = get<0>(odesol) ;
	const vector<double>& y = get<1>(odesol) ;

	ofstream out ;
	out.open (filename, std::ios::out) ;

	if ( !out.is_open() ) {
		throw runtime_error("I/O Error: " + filename) ;
	}

	streamsize prec = out.precision() ;
	typedef std::vector<double>::size_type size ;
	for (size i = 0 ; i != t.size() ; ++i) {
		out << scientific << setprecision(15) << t[i] << "\t" 
		     << y[i] << setprecision(prec) << endl ;
	}

	out.close() ;
}
