/**
* \file ComputeGaussKronrodNodesWeights.h
* The functions contained in this file calculate the Gauss-Kronrod nodes and weights
* using the Laurie/Gautschi method.
*/

#ifndef EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H
#define EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H

#include "LaurieGautschi.h" // stable, slow and requires c++11 tgamma
#include "Piessens.h" // stable for most rules and precisions, faster than LaurieGautschi
#include "Monegato.h" // fast but unstable at the moment for rules above 80

#endif //EIGEN_COMPUTE_GAUSS_KRONROD_NODES_WEIGHTS_H
