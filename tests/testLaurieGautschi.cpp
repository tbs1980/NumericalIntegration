#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

int compare_codes(void)
{
    typedef double RealType;
    //typedef mpfr::mpreal RealType;
    //RealType::set_default_prec(128);
    //RealType::set_default_prec(256);

    std::cout<<std::endl<<"MS Laurie Gautschi"<<std::endl;
    Eigen::Array<RealType,Dynamic,2> ans;
    int n = 10;
    ans = Kronrod::multiPrecisionKronrod<RealType>(n);

    std::cout<<std::fixed;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout<<std::setprecision(15)<<ans(i,0)<<"\t"<<ans(i,1)<<std::endl;
    }

    std::cout<<"\nSTB Laurie Gautschi"<<std::endl;

    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N=10;
    VectorType x=VectorType::Zero(2*N+1);
    VectorType w=VectorType::Zero(2*N+1);

    LaurieGautschiPolicy::mpkonrad(N,x,w);

    for(IndexType i=0;i<x.rows();++i)
    {
        std::cout<<std::setprecision(15)<<x(i)<<"\t"<<w(i)<<std::endl;
    }

    Eigen::Array<RealType, Dynamic, 1> xGKPosAndNeg = Eigen::Array<RealType, Dynamic, 1>::Zero(x.rows());
    Eigen::Array<RealType, Dynamic, 1> wGKPosAndNeg = Eigen::Array<RealType, Dynamic, 1>::Zero(x.rows());;
    Eigen::Array<RealType, Dynamic, 1> xGK;
    Eigen::Array<RealType, Dynamic, 1> wGK;
    Eigen::Array<RealType, Dynamic, 1> wG;

    Kronrod::kronrod(n, xGK,  wGK, wG);

    for(IndexType i=0; i<xGK.rows(); ++i)
    {
        xGKPosAndNeg(i) = -xGK(i);
        wGKPosAndNeg(i) = wGK(i);
    }

    for(IndexType i=0; i<xGK.rows(); ++i)
    {
        xGKPosAndNeg(xGKPosAndNeg.rows()-1-i) = xGK(i);
        wGKPosAndNeg(wGKPosAndNeg.rows()-1-i) = wGK(i);
    }

    std::cout << "\nMS Piessens" << std::endl;
    for(IndexType i=0;i<xGKPosAndNeg.rows();++i)
    {
        std::cout<<std::setprecision(15)<<xGKPosAndNeg(i)<<"\t"<<wGKPosAndNeg(i)<<std::endl;
    }

    std::cout<<"\nSolution Differences: MS - STB Laurie Gautschi"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout<<std::setprecision(15)<<ans(i,0) - x(i)<<"\t"<<ans(i,1) - w(i)<<std::endl;
    }

    std::cout<<"\nSolution Differences: MS Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout<<std::setprecision(15)<<ans(i,0) - xGKPosAndNeg(i)<<"\t"<<ans(i,1) - wGKPosAndNeg(i)<<std::endl;
    }

    std::cout<<"\nSolution Differences: STB Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<x.rows();++i)
    {
        std::cout<<std::setprecision(15)<< x(i)-xGKPosAndNeg(i)<<"\t"<<w(i)-wGKPosAndNeg(i)<<std::endl;
    }

    std::cout<<std::endl;

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=compare_codes();

    return ret;
}
