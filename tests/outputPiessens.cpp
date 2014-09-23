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

    const IndexType N=30;
    const int outputIntegers = 34;

    Eigen::Array<RealType,Dynamic,2> ans;
    ans = Kronrod::multiPrecisionKronrod<RealType>(N);

    std::cout<<std::fixed<<std::endl<<"Rule: "<<2*N+1<<std::endl;

    Eigen::Array<RealType, Dynamic, 1> xGK;
    Eigen::Array<RealType, Dynamic, 1> wGK;
    Eigen::Array<RealType, Dynamic, 1> wG;

    Kronrod::kronrod(N, xGK,  wGK, wG);

    for(IndexType i=0;i<xGK.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << xGK(i) << "," << std::endl;
    }

    std::cout<<std::endl;

    for(IndexType i=0;i<wGK.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << wGK(i) << "," << std::endl;
    }

    std::cout << std::endl;

    for(IndexType i=0;i<wG.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << wG(i) << "," << std::endl;
    }

    std::cout<<std::endl;

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=compare_codes();

    return ret;
}
    const IndexType N=7;
    const int outputIntegers = 33;

    Eigen::Array<RealType,Dynamic,2> ans;
    ans = Kronrod::multiPrecisionKronrod<RealType>(N);

    std::cout<<std::fixed<<std::endl<<"Rule: "<<2*N+1<<std::endl;

    Eigen::Array<RealType, Dynamic, 1> xGK;
    Eigen::Array<RealType, Dynamic, 1> wGK;
    Eigen::Array<RealType, Dynamic, 1> wG;

    Kronrod::kronrod(N, xGK,  wGK, wG);

    for(IndexType i=0;i<xGK.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << xGK(i) << "," << std::endl;
    }

    std::cout<<std::endl;

    for(IndexType i=0;i<wGK.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << wGK(i) << "," << std::endl;
    }

    std::cout << std::endl;

    for(IndexType i=0;i<wG.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << "    " << wG(i) << "," << std::endl;
    }

    std::cout<<std::endl;

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=compare_codes();

    return ret;
}
