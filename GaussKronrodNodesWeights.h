/**
* \file GaussKronrodNodesWeights.h
* The functions contained in this file calculate the Gauss-Kronrod nodes and weights
* using the Laurie/Gautschi method.
*/

#ifndef NI_KRONRODLAURIEGAUTSCHI_H
#define NI_KRONRODLAURIEGAUTSCHI_H

#include "LaurieGautschi.h" // stable, slow and requires c++11 tgamma
#include "Piessens.h" // stable for most rules and precisions, faster than LaurieGautschi
#include "Monegato.h" // fast but unstable at the moment for rules above 80

#endif //NI_KRONRODLAURIEGAUTSCHI_H
