/**
 * \file testNodesAndWeightsByCompare.cpp
 * This file is a unit test for ComputeGaussKronrodNodesWeights.h.
 * The test is used to compare two different approaches to computing Gauss-Kronrod node and
 * weight values for a single ruleset for agreement at the specified level of precision.
 */

#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

int compare_codes_unified_interface(const unsigned int N=10)
{
    // typedef float Scalar;
    // typedef double Scalar;
    // typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(320); // 128, 320, 384, 448 gives an error Newton-Raphson iterative abscissae solver failed.
    // 256,288,320,352,384,416,448


    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGLaurieGautschi;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGLaurieGautschi;

    typedef Eigen::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef Eigen::Piessens<Scalar> PiessensPolicy;
    typedef Eigen::Monegato<Scalar> MonegatoPolicy;

    const int outputIntegers = 20; //beyond 67 integers methods disagree for precision=256.

    LaurieGautschiPolicy::computeAbscissaeAndWeights(N, xGKLaurieGautschi, wGKLaurieGautschi, xGLaurieGautschi, wGLaurieGautschi);

    std::ofstream fout;

    std::string fileLocation = "test/testOutput/";
    std::string fileName = "LaurieGautschi320.dat";
    std::string fileNameAndLocation = fileLocation + fileName;
    fout.open(fileNameAndLocation);

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << std::endl << "Kronrod Nodes" << std::endl;
    fout << std::fixed;
    
    for(DenseIndex i = 0; i < xGKLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKLaurieGautschi(i) << ",\n";
    }

    fout << std::endl << "Kronrod Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGKLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKLaurieGautschi(i) << ",\n";
    }

    fout << std::endl << "Gauss Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGLaurieGautschi.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGLaurieGautschi(i) << ",\n";
    }

    fout.close();

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGPiessens;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGPiessens;

    PiessensPolicy::computeAbscissaeAndWeights(N,xGKPiessens,wGKPiessens,xGPiessens,wGPiessens);

    fileName = "Piessens320.dat";
    fileNameAndLocation = fileLocation + fileName;
    fout.open(fileNameAndLocation);

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << std::endl << "Kronrod Nodes" << std::endl;
    fout << std::fixed;
    
    for(DenseIndex i = 0; i < xGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKPiessens(i) << ",\n";
    }

    fout << std::endl <<"Kronrod Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGKPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKPiessens(i) << ",\n";
    }

    fout << std::endl << "Gauss Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGPiessens(i) << ",\n";
    }

    fout.close();

    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGKMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> xGMonegato;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> wGMonegato;

    MonegatoPolicy::computeAbscissaeAndWeights(N,xGKMonegato,wGKMonegato,xGMonegato,wGMonegato);

    fileName = "Monegato320.dat";
    fileNameAndLocation = fileLocation + fileName;
    fout.open(fileNameAndLocation);

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << std::endl << "Kronrod Nodes" << std::endl;
    fout << std::fixed;
    
    for(DenseIndex i = 0; i < xGKMonegato.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGKMonegato(i) << ",\n";
    }

    fout << std::endl << "Kronrod Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGKMonegato.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGKMonegato(i) << ",\n";
    }

    fout << std::endl << "Gauss Weights" << std::endl;
    
    for(DenseIndex i = 0; i < wGMonegato.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGMonegato(i) << ",\n";
    }

    fout.close();

    return EXIT_SUCCESS;
}


int main(int argc, char** argv)
{
    int m = (argc > 1) ? atoi(argv[1]) : 10;	// Legendre degree
    
    int ret = EXIT_SUCCESS;
    
    ret += compare_codes_unified_interface(m);
    
    return ret;
}
