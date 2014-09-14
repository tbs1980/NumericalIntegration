#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

int main(void)
{
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);

    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N=10;
    VectorType x=VectorType::Zero(2*N+1);
    VectorType w=VectorType::Zero(2*N+1);

    LaurieGautschiPolicy::mpkonrad(N,x,w);

    std::cout<<std::fixed;
    for(IndexType i=0;i<x.rows();++i)
    {
        std::cout<<std::setprecision(25)<<x(i)<<"\t"<<w(i)<<std::endl;
    }

    return 0;
}
