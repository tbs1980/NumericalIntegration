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
 * \defgroup NumericalIntegration_Module
 * \brief This module provides an adaptive quadrature method of numerical integration.
 *
 * This module provides an adaptive quadrature method of numerical integration of the style
 * implemented in the QUADPACK library while offering functionality to calculate nodes/weights
 * for Gauss-Kronrod integration, unit tests, and support for multiprecision using mpreal
 * precision type.
 *
 * To use this module, add
 * \code
 * #include <unsupported/Eigen/NumericalIntegration>
 * \endcode
 * at the start of your source file.
 */
}

#include <iomanip>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MPRealSupport>

#include "ComputeGaussKronrodNodesWeights.h"
#include "GaussKronrodNodesWeights.h"
#include "Integrator.h"

#endif // EIGEN_NUMERICAL_INTEGRATION_H
