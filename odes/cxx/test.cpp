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
using ode::EulerRK2 ;


void write (std::tuple < std::vector<double>, std::vector<double> > ) ;
void write (std::tuple < std::vector<double>, std::vector<double> >,
	    const std::string& ) ;

int main() {

	const int N = 255;
	const double yi = 1.0 ;
	const double ti = 0.0, tf = 5.0 ;
	tuple < vector<double>, vector<double> > odesol_Euler, odesol_RK2 ;

	// odefun
	auto f = [](const double& t, const double& y) {
		const double k = 1.0 ;
		return (-k * y) ;
       	} ;

	// solves the ODE using the specified method
	odesol_Euler = Euler    (ti, tf, yi, N, f) ;
	odesol_RK2   = EulerRK2 (ti, tf, yi, N, f) ;
	
	// exports numerical solutions
	string filename = "output/Euler.dat" ;
	write (odesol_Euler, filename) ;

	filename = "output/EuRK2.dat" ;
	write (odesol_RK2, filename) ;

	write (odesol_RK2) ;
	return 0 ;
}


void write (std::tuple < std::vector<double>, std::vector<double> > odesol)
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


void write (std::tuple < std::vector<double>, std::vector<double> > odesol,
	    const std::string& filename)
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
