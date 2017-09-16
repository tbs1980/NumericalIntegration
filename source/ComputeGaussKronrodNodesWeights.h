/**
* \file ComputeGaussKronrodNodesWeights.h
* The functions contained in this file calculate the Gauss-Kronrod nodes and weights
* using the Laurie/Gautschi method.
*/

#ifndef EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H
#define EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H

#include "LaurieGautschi.h" // Stable, slow and requires c++11 tgamma and the latest mpreal.h.
#include "Piessens.h" 		// Stable for most rules and precisions, faster than LaurieGautschi.
#include "Monegato.h" 		// Fastest but unstable for rules above 80.

#endif //EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H
