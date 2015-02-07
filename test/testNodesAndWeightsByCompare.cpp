/**
 * \file testNodesAndWeightsByCompare.cpp
 * This file is a unit test for GaussKronrodNodesWeights.h.
 * The test is used to compare two different approaches to computing Gauss-Kronrod node and
 * weight values for a single ruleset for agreement at the specified level of precision.
 */

#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

int compare_codes_unified_interface(void)
{
    //typedef float Scalar;
    typedef double Scalar;
    //typedef long double Scalar;
    //typedef mpfr::mpreal Scalar;
    //Scalar::set_default_prec(320); //128, 320, 384, 448 gives an error Newton-Raphson iterative abscissae solver failed.
    //256,288,320,352,384,416,448
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGLaurieGautschi;

    typedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef Kronrod::Piessens<Scalar> PiessensPolicy;
    typedef Kronrod::Monegato<Scalar> MonegatoPolicy;

    typedef LaurieGautschiPolicy::IndexType IndexType;

    const unsigned int N = 100;
    const int outputIntegers = 80; //beyond 67 integers methods disagree for precision=256.

    LaurieGautschiPolicy::computeAbscissaeAndWeights(N,xGKLaurieGautschi,wGKLaurieGautschi,xGLaurieGautschi,wGLaurieGautschi);

    std::ofstream fout;
    fout.open("LaurieGautschi320.dat");

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << "\nKronrod Nodes\n";
    fout << std::fixed;
    for(IndexType i = 0; i < xGKLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKLaurieGautschi(i) << ",\n";
    }

    fout << "\nKronrod Weights\n";
    for(IndexType i = 0; i < wGKLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKLaurieGautschi(i) << ",\n";
    }

    fout << "\nGauss Weights\n";
    for(IndexType i = 0; i < wGLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGLaurieGautschi(i) << ",\n";
    }

    fout.close();

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGPiessens;

    PiessensPolicy::computeAbscissaeAndWeights(N,xGKPiessens,wGKPiessens,xGPiessens,wGPiessens);

    fout.open("Piessens320.dat");

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << "\nKronrod Nodes\n";
    fout << std::fixed;
    for(IndexType i = 0; i < xGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKPiessens(i) << ",\n";
    }

    fout << "\nKronrod Weights\n";
    for(IndexType i = 0; i < wGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKPiessens(i) << ",\n";
    }

    fout << "\nGauss Weights\n";
    for(IndexType i = 0; i < wGPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGPiessens(i) << ",\n";
    }

    fout.close();

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGMonegato;

    MonegatoPolicy::computeAbscissaeAndWeights(N,xGKMonegato,wGKMonegato,xGMonegato,wGMonegato);

    fout.open("Monegato320.dat");

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << "\nKronrod Nodes\n";
    fout << std::fixed;
    for(IndexType i = 0; i < xGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKPiessens(i) << ",\n";
    }

    fout << "\nKronrod Weights\n";
    for(IndexType i = 0; i < wGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKPiessens(i) << ",\n";
    }

    fout << "\nGauss Weights\n";
    for(IndexType i = 0; i < wGPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGPiessens(i) << ",\n";
    }

    fout.close();

    return EXIT_SUCCESS;
}


int main(void)
{
    int ret = 0;
    ret += compare_codes_unified_interface();
    return ret;
}
