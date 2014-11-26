#include <NIHeaders.h>

#include <iostream>
#include <fstream>
#include <iomanip>

int compare_codes(void)
{
    std::ofstream fout;
    fout.open("SingleRuleKronrodNodesAndWeights.txt");

    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);

    typedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N = 100;
    const int outputIntegers = 256;

    VectorType xGK = VectorType::Zero(2*N+1);
    VectorType wGK = VectorType::Zero(2*N+1);
    VectorType xG = VectorType::Zero(N);
    VectorType wG = VectorType::Zero(N);

    LaurieGautschiPolicy::mpkronrod(N,xGK,wGK);
    LaurieGautschiPolicy::mpgauss(N,xG,wG);

    fout << "Kronrod Nodes and Weights for N = " << N << std::endl;
    fout << "\nLaurie Gautschi Calculations\n";
    fout << "\nKronrod Nodes\n";
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xGK(i) << ",\n";
    }

    fout << "\nKronrod Weights\n";
    for(IndexType i = 0; i < wGK.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGK(i) << ",\n";
    }

    fout << "\nGauss Nodes\n";
    for(IndexType i = 0; i < xG.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << xG(i) << ",\n";
    }
    
    fout << "\nGauss Weights\n";
    for(IndexType i = 0; i < xG.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wG(i) << ",\n";
    }

    Eigen::Array<Scalar, Dynamic, 1> xGKPiessens;
    Eigen::Array<Scalar, Dynamic, 1> wGKPiessens;
    Eigen::Array<Scalar, Dynamic, 1> wGPiessens;

    Kronrod::kronrod(N, xGKPiessens,  wGKPiessens, wGPiessens);

    for(IndexType i = 0; i < xGKPiessens.rows(); ++i)
    {
        xGK(i) = -xGKPiessens(i);
        wGK(i) = wGKPiessens(i);
    }

    for(IndexType i=0; i<xGKPiessens.rows(); ++i)
    {
        xGK(xGK.rows()-1-i) = xGKPiessens(i);
        wGK(wGK.rows()-1-i) = wGKPiessens(i);
    }

    fout << "\nPiessens Calculations\n";
    fout << "\nKronrod Nodes" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
       fout << std::setprecision(outputIntegers) << xGK(i) << ",\n";
    }
    
    fout << "\nKronrod Weights" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGK(i) << ",\n";
    }

    fout << "\n(Piessens' method does not calculate the Gauss Nodes)\n";
    fout << "\nGauss Weights" << std::endl;
    for(IndexType i = 0; i < wGPiessens.rows(); ++i)
    {
        fout << std::setprecision(outputIntegers) << wGPiessens(i) << ",\n";
    }

    fout.close();

    std::cout << "\nNode Calculations successfully written to \"SingleRuleKronrodNodesAndWeights.txt\"." << std::endl;

    return EXIT_SUCCESS;
}


int main(void)
{
    int ret = 0;
    ret += compare_codes();

    return ret;
}
