#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

int compare_codes(void)
{
    std::cout<<"MS"<<std::endl;
    Eigen::Array<double,Dynamic,2> ans;
    int n = 10;
    ans = Kronrod::multiPrecisionKronrod(n);

    std::cout<<std::fixed;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout<<std::setprecision(15)<<ans(i,0)<<"\t"<<ans(i,1)<<std::endl;
    }

    std::cout<<"\nSTB"<<std::endl;

    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    //RealType::set_default_prec(256);

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

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=compare_codes();

    return ret;
}
