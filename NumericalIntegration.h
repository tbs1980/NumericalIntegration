// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Development work based on work from QUADPACK, Robert Piessens, et al,
// and original work by Dirk Laurie, Walter Gautschi, with support by work
// from John Burkardt.
//
// Multiprecision templating by Pavel Holoborodko and
// code porting, multiprecision templating, and unit tests created by
// Sreekumar Thaithara Balan, Mark Sauder, and Matt Beall 2014.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_NUMERICAL_INTEGRATION_H
#define EIGEN_NUMERICAL_INTEGRATION_H

namespace Eigen
{
/**
 * \defgroup Numerical_Integration_Module Quadrature and Nodes_Weights module
 *
 * This module provides an adaptive quadrature method of numerical integration of the style
 * implemented in the QUADPACK library while offering functionality to calculate nodes/weights
 * for Gauss-Kronrod integration, unit tests, and support for multiprecision using mpreal
 * precision type.
 */

    //include <unsupported/Eigen/NumericalIntegration>
}

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MPRealSupport>

//http://www.wolframalpha.com/input/?i=pi+to+500+decimal+places
//http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
#define NI_M_PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194913

// \TODO Resolve the following overload requirements within mpfr
// Need to switch between gamma() from mpfrc++ and tgamma() from c++
template<typename T>
T Gamma(T x)
{
	return gamma(x);
}

double Gamma(double x)
{
	return tgamma(x);
}

float Gamma(float x)
{
	return tgamma(x);
}

long double Gamma(long double x)
{
	return tgamma(x);
}

#include "KronrodLaurieGautschi.h"
#include "KronrodPiessens.h"

#include "QuadratureKronrod.h"
#include "Integrator.h"

#endif // EIGEN_NUMERICAL_INTEGRATION_H
