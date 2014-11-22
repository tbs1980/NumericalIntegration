#include <NIHeaders.h>
#include <iostream>
#include <iomanip>

int compare_codes(void)
{
    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);

    typedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N = 10;
    const int outputIntegers = 33;

    VectorType xGK = VectorType::Zero(2*N+1);
    VectorType wGK = VectorType::Zero(2*N+1);
    VectorType xG = VectorType::Zero(N);
    VectorType wG = VectorType::Zero(N);

    LaurieGautschiPolicy::mpkronrod(N,xGK,wGK);
    LaurieGautschiPolicy::mpgauss(N,xG,wG);

    std::cout << "\nLaurie Gautschi Calculations" << std::endl;
    std::cout << "\nKronrod Nodes" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << xGK(i) << std::endl;
    }

    std::cout << "\nKronrod Weights" << std::endl;
    for(IndexType i = 0; i < wGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << wGK(i) << std::endl;
    }

    std::cout << "\nGauss Nodes" << std::endl;
    for(IndexType i = 0; i < xG.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << xG(i) << std::endl;
    }
    
    std::cout << "\nGauss Weights" << std::endl;
    for(IndexType i = 0; i < xG.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << wG(i) << std::endl;
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

    std::cout << "\nPiessens Calculations" << std::endl;
    std::cout << "\nKronrod Nodes" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << xGK(i) << std::endl;
    }
    
    std::cout << "\nKronrod Weights" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << wGK(i) << std::endl;
    }

    std::cout << "\n(Piessens' method does not calculate the Gauss Nodes)" << std::endl;
    std::cout << "\nGauss Weights" << std::endl;
    for(IndexType i = wGPiessens.rows() - 1; i >= 0; --i)
    {
        std::cout << std::setprecision(outputIntegers) << wGPiessens(i) << std::endl;
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}


int main(void)
{
    int ret = 0;
    ret += compare_codes();

    return ret;
}
