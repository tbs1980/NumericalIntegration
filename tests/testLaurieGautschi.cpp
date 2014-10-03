#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

int compare_codes(void)
{
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    //RealType::set_default_prec(128);
    RealType::set_default_prec(256);

    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N = 10;
    const int outputIntegers = 33;

    Eigen::Array<RealType,Dynamic,2> xwGK;
    Eigen::Array<RealType,Dynamic,2> xwG;

    xwGK = Kronrod::multiPrecisionKronrod<RealType>(N);
    xwG = Kronrod::multiPrecisionGauss<RealType>(N);

    std::cout<<std::fixed;
    std::cout<<std::endl<<"MS Laurie Gautschi"<<std::endl;
    for(int i = 0; i < xwGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers);
        std::cout << xwGK.col(0)(i) << "\t" << xwGK.col(1)(i) << std::endl;
    }

    std::cout<<std::endl;
    for(int i = 0; i < xwG.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers);
        std::cout << xwG.col(0)(i) << "\t" << xwG.col(1)(i) << std::endl;
    }

    VectorType xGK =VectorType::Zero(2*N+1);
    VectorType wGK =VectorType::Zero(2*N+1);
    LaurieGautschiPolicy::mpkonrad(N,xGK,wGK);

    std::cout<<"\nSTB Laurie Gautschi"<<std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers);
        std::cout << xGK(i) << "\t" << wGK(i) << std::endl;
    }

    VectorType xG = VectorType::Zero(N);
    VectorType wG = VectorType::Zero(N);

    LaurieGautschiPolicy::mpgauss(N,xG,wG);

    std::cout<<std::endl;
    for(IndexType i = 0; i < xG.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers);
        std::cout << xG(i) << "\t" << wG(i) << std::endl;
    }

    Eigen::Array<RealType, Dynamic, 1> xGKPiessens;
    Eigen::Array<RealType, Dynamic, 1> wGKPiessens;
    Eigen::Array<RealType, Dynamic, 1> wGPiessens;

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

    std::cout << "\nMS Piessens" << std::endl;
    for(IndexType i = 0; i < xGK.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers);
        std::cout << xGK(i) << "\t" << wGK(i) << std::endl;
    }

    std::cout << std::endl;

    for(IndexType i = 0; i < wGPiessens.rows(); ++i)
    {
        std::cout << std::setprecision(outputIntegers) << "\t\t\t\t\t";
        std::cout << wGPiessens(i) << std::endl;
    }

    for(IndexType i = wGPiessens.rows() - 1; i >= 0; --i)
    {
        std::cout << std::setprecision(outputIntegers) << "\t\t\t\t\t";
        std::cout << wGPiessens(i) << std::endl;
    }

    std::cout<<std::endl;

/*------------------------Output of Differences Between Approaches-------------------------------//

    std::cout<<"\nSolution Differences: MS - STB Laurie Gautschi"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout << std::setprecision(15) << ans(i,0) - x(i)
          << "\t" << ans(i,1) - w(i) << std::endl;
    }

    std::cout<<"\nSolution Differences: MS Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout << std::setprecision(15) << ans(i,0) - xGKPosAndNeg(i)
          << "\t" << ans(i,1) - wGKPosAndNeg(i) << std::endl;
    }

    std::cout<<"\nSolution Differences: STB Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<x.rows();++i)
    {
        std::cout  << std::setprecision(15) << x(i) - xGKPosAndNeg(i)
           << "\t" << w(i) - wGKPosAndNeg(i) << std::endl;
    }

//-----------------------End Output of Differences Between Approaches----------------------------*/

    return EXIT_SUCCESS;
}


int main(void)
{
    int ret = 0;
    ret += compare_codes();

    return ret;
}
