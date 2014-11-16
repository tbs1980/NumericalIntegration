// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Development work based on work from QUADPACK, Robert Piessens, et al,
// and original work by Dirk Laurie, Walter Gautschi, with support by work
// from John Burkardt. Code porting, multiprecision templating, and unit 
// tests created by Pavel Holoborodko, Sreekumar Thaithara Balan, Mark Sauder,
// and Matt Beall 2014
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

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MPRealSupport>

#include "nodes_weights/kronrodLaurieGautschi.h"
#include "nodes_weights/kronrodPiessens.h"

#include "quadrature/QuadratureKronrod.h"
#include "quadrature/Integrator.h"

#endif // EIGEN_NUMERICAL_INTEGRATION_H

