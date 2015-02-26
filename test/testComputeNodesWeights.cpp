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

int test_nodes_weights_difference(const unsigned int N)
{
    //typedef float Scalar;
    //typedef double Scalar;
    typedef long double Scalar;
    //typedef mpfr::mpreal Scalar;
    //Scalar::set_default_prec(50);

    typedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef Kronrod::Piessens<Scalar> PiessensPolicy;
    typedef Kronrod::Monegato<Scalar> MonegatoPolicy;

    typedef LaurieGautschiPolicy::IndexType IndexType;

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

    using std::abs;
    for(IndexType i = 0; i < xGKLaurieGautschi.rows(); ++i)
    {
        assert(abs(xGKLaurieGautschi(i) - xGKPiessens(i)) > 1e-200);
        assert(abs(xGKLaurieGautschi(i) - xGKMonegato(i)) > 1e-200);
        assert(abs(xGKPiessens(i) - xGKMonegato(i)) > 1e-200);
    }

    for(IndexType i = 0; i < wGKLaurieGautschi.rows(); ++i)
    {
        assert((wGKLaurieGautschi(i) - wGKPiessens(i)) > 1e-50);
        assert((wGKLaurieGautschi(i) - wGKMonegato(i)) > 1e-50);
        assert((wGKPiessens(i) - wGKMonegato(i)) > 1e-50);
    }

    for(IndexType i = 0; i < xGLaurieGautschi.rows(); ++i)
    {
        assert((xGLaurieGautschi(i) - xGPiessens(i)) > 1e-50);
        assert((xGLaurieGautschi(i) - xGMonegato(i)) > 1e-50);
        assert((xGPiessens(i) - xGMonegato(i)) > 1e-50);
    }

    for(IndexType i = 0; i < wGLaurieGautschi.rows(); ++i)
    {
        assert((wGLaurieGautschi(i) - wGPiessens(i)) > 1e-50);
        assert((wGLaurieGautschi(i) - wGMonegato(i)) > 1e-50);
        assert((wGPiessens(i) - wGMonegato(i)) > 1e-50);
    }

    return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
    size_t m = (argc > 1) ? atoi(argv[1]) : 10;	// Legendre degree

    int ret = EXIT_SUCCESS;
    ret += test_nodes_weights_difference(m);
    return ret;
}
