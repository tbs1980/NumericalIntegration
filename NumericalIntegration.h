// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Development work based on work from QUADPACK, Robert Piessens, et al,
// and original work by Dirk Laurie, Walter Gautschi, with support by work
// from John Burkardt.
//
// Multiprecision templating by Pavel Holoborodko.
// Code porting and unit tests created by Sreekumar Thaithara Balan,
// Mark Sauder, and Matt Beall 2014
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
  *
  * \code
  * #include <unsupported/Eigen/NumericalIntegration>
  * \endcode
  */
}
//#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MPRealSupport>

// @TODO resolve the  following overload requirements

template<typename T>
T Pow(T base, T exponent)
{
    return pow(base,exponent);
}

float Pow(float base, float exponent)
{
    return std::pow(base,exponent);
}

double Pow(double base, double exponent)
{
    return std::pow(base,exponent);
}

long double Pow(long double base, long double exponent)
{
    return std::pow(base,exponent);
}

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

template<typename T>
T Sin(T arg)
{
    return sin(arg);
}

float Sin(float arg)
{
    return std::sin(arg);
}

double Sin(double arg)
{
    return std::sin(arg);
}

long double Sin(long double arg)
{
    return std::sin(arg);
}

template<typename T>
T Abs(T arg)
{
    return abs(arg);
}

float Abs(float arg)
{
    return std::abs(arg);
}

double Abs(double arg)
{
    return std::abs(arg);
}

long double Abs(long double arg)
{
    return std::abs(arg);
}


#include "nodes_weights/kronrodLaurieGautschi.h"
#include "nodes_weights/kronrodPiessens.h"
#include "nodes_weights/kronrodPiessensClass.h"

#include "quadrature/QuadratureKronrod.h"
#include "quadrature/Integrator.h"

#endif // EIGEN_NUMERICAL_INTEGRATION_H
