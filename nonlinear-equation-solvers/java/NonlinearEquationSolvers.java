/*
 * Computational Solutions                               May 01, 2022
 * Prof. M. Diaz-Maldonado
 *
 *
 * Synopsis:
 * Recursive implementation of (some) nonlinear equation solvers.
 *
 *
 * Copyright (c) 2022 Misael Diaz-Maldonado
 * This file is released under the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 *
 * References:
 * [0] JJ McConnell, Analysis of Algorithms, 2nd edition.
 * [1] A Gilat and V Subramanian, Numerical Methods for Engineers and
 *     Scientists, 3rd edition.
 * [2] BR Munson et al., Fundamentals of Fluid Mechanics, 8th edition.
 *
 */

import static java.lang.Math.sqrt;
import static java.lang.Math.log10;

public class NonlinearEquationSolvers
{

  interface Objf
  {
    public double f (double x);			// nonlinear function, f(x)
  }


  // defines the nonlinear objective function f(x) as a lambda
  static Objf objf = (double x) -> {
    return (1.0/sqrt(x) + 2.0 * log10(0.024490 / 3.7 + 2.51 / (7812525.1 * sqrt(x) ) ) );
  };


  public static void main (String [] args)	// solves for nonlinear function f(x):
  {
    Bisection();				// applies recursive Bisection Method
    RegulaFalsi();				// applies recursive Regula Falsi Method
  }


  // double bisect (double l, double u, Objf objf)
  //
  // Synopsis:
  // Recursive implementation of the Bisection Method.
  //
  // Inputs
  // l		lower bracketing limit
  // u 		upper bracketing limit
  // objf	nonlinear function f(x)
  //
  // Output
  // x		returns root of the nonlinear function f(x)


  public static double bisect (double l, double u, Objf objf)
  {
    hasRoot(l, u, objf);			// complains if [l, u] bracket no root
    double tol = 1.0e-12;			// tolerance
    double m = (l + u) / 2;			// middle value
    double f;					// f(x_m)

    f = (objf.f(m) < 0)? -objf.f(m) : objf.f(m);// gets the absolute value of f(x_m)

    if (f < tol)				// if f(x_m) meets convergence criterion:
    {
      return m;					// uses direct solution
    }
    else					// divides otherwise:
    {
      if (objf.f(m) * objf.f(u) < 0.0)
      {
	l = m;					// divides by removing first interval
      }
      else
      {
	u = m;					// divides by removing second interval
      }
    }

    return bisect(l, u, objf);			// recurses
  }


  // double regfal (double l, double u, Objf objf)
  //
  // Synopsis:
  // Recursive implementation of the Regula Falsi (or False Position) Method.
  //
  // Inputs
  // l		lower bracketing limit
  // u 		upper bracketing limit
  // objf	nonlinear function f(x)
  //
  // Output
  // x		returns root of the nonlinear function f(x)


  public static double regfal (double l, double u, Objf objf)
  {
    hasRoot(l, u, objf); 			// complains if [l, u] brackets no root
    double tol = 1.0e-12;			// tolerance
    double fl = objf.f(l), fu = objf.f(u);	// gets f(x_l) and f(x_u)
    double m = (l * fu - u * fl) / (fu - fl); 	// interpolates to get intermediate value 
    double f;					// f(x_m)

    f = (objf.f(m) < 0)? -objf.f(m) : objf.f(m);// gets the absolute value of f(x_m)

    if (f < tol)				// if f(x_m) meets convergence criterion:
    {
      return m;					// uses direct solution
    }
    else					// divides otherwise:
    {
      if (objf.f(m) * objf.f(u) < 0.0)
      {
	l = m;					// divides by removing first interval
      }
      else
      {
	u = m;					// divides by removing second interval
      }
    }

    return regfal(l, u, objf);			// recurses
  }


  private static void hasRoot (double l, double u, Objf objf)
  {
    // complains if the bracketing interval [l, u] does not bracket a root
    if ( objf.f(l) * objf.f(u) > 0.0 )
    {
      String e = "solver(): expects a root in the interval [l, u]";
      throw new IllegalArgumentException(e);
    }
  }


  private static void Bisection ()		// solves f(x) with Bisection
  {
    double l = 2.0e-2, u = 6.0e-2;		// defines the bracketing interval [l, u]
    double r = bisect (l, u, objf);		// solves for the root of f(x) recursively
    System.out.println("\nBisection Method:\n");
    System.out.printf("root: %.12f\n", r);
    System.out.printf("f(x): %.12e\n", objf.f(r) );
  }


  private static void RegulaFalsi ()		// solves f(x) with Regula Falsi
  {
    double l = 2.0e-2, u = 6.0e-2;		// defines the bracketing interval [l, u]
    double r = regfal (l, u, objf);		// solves for the root of f(x) recursively
    System.out.println("\nRegula Falsi Method:\n");
    System.out.printf("root: %.12f\n", r);
    System.out.printf("f(x): %.12e\n", objf.f(r) );
  }
}
