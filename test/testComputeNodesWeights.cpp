/**
 * \file testComputeNodesWeights.cpp
 * This file is a unit test for ComputeGaussKronrodNodesWeights.h.
 * The test is used to compare two different approaches to computing Gauss-Kronrod node and
 * weight values for a single ruleset for agreement at the specified level of precision.
 */

#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

using namespace Eigen;

int test_nodes_weights_difference(const unsigned int N)
{
    //typedef float Scalar;
    //typedef double Scalar;
    // typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(50);

    typedef Eigen::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef Eigen::Piessens<Scalar> PiessensPolicy;
    typedef Eigen::Monegato<Scalar> MonegatoPolicy;

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGLaurieGautschi;

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGPiessens;

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGMonegato;

    LaurieGautschiPolicy::computeAbscissaeAndWeights(N,xGKLaurieGautschi,wGKLaurieGautschi,xGLaurieGautschi,wGLaurieGautschi);
    PiessensPolicy::computeAbscissaeAndWeights(N,xGKPiessens,wGKPiessens,xGPiessens,wGPiessens);
    MonegatoPolicy::computeAbscissaeAndWeights(N,xGKMonegato,wGKMonegato,xGMonegato,wGMonegato);

    double epsilon = 1e-16;

    for(Index i = 0; i < xGKLaurieGautschi.rows(); ++i)
    {
        if ((abs(xGKLaurieGautschi(i) - xGKPiessens(i)) > epsilon)
           || (abs(xGKLaurieGautschi(i) - xGKMonegato(i)) > epsilon)
           || (abs(xGKPiessens(i) - xGKMonegato(i)) > epsilon))
        {
            std::cout << "Failed xGK " << i << std::endl;
            return EXIT_FAILURE;
        }
    }

    for(Index i = 0; i < wGKLaurieGautschi.rows(); ++i)
    {
        if ((abs(wGKLaurieGautschi(i) - wGKPiessens(i)) > epsilon)
           || (abs(wGKLaurieGautschi(i) - wGKMonegato(i)) > epsilon)
           || (abs(wGKPiessens(i) - wGKMonegato(i)) > epsilon))
        {
            std::cout << "Failed wGK " << i << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    // Note: Piessens and Monegato methods do not calculate the gauss nodes independently.
    //       Piessens and Monegato gauss nodes are taken directly from GaussKronrod nodes,
    //       and agreement in the GaussKronrod node test implies agreement for Gauss nodes.

    for(Index i = 0; i < wGLaurieGautschi.rows(); ++i)
    {
        if ((abs(wGLaurieGautschi(i) - wGPiessens(i)) > epsilon)
           || (abs(wGLaurieGautschi(i) - wGMonegato(i)) > epsilon)
           || (abs(wGPiessens(i) - wGMonegato(i)) > epsilon))
        {
            std::cout << "Failed wG " << i << std::endl;
            return EXIT_FAILURE;
        }
    }

    std::cout << "\n\tTests Succeeded!\n" << std::endl;
    return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
    int m = (argc > 1) ? atoi(argv[1]) : 10;	// Legendre degree

    int ret = EXIT_SUCCESS;
    
    ret += test_nodes_weights_difference(m);
    return ret;
}
