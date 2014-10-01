#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>
#include "../quadrature/QuadratureKronrod_STB.hpp"

int test_values()
{
    typedef float RealType;
    //typedef double RealType;
    //typedef mpfr::mpreal RealType;
    //RealType::set_default_prec(256);
    typedef Eigen::QuadratureKronrod<RealType> QuadratureKronrodValuesType;

    //QuadratureKronrodValuesType qgk;
    QuadratureKronrodValuesType::ComputeNodesAndWeights();

    std::cout<<"\nGaussKronrod15 \n"<<std::endl;
    std::cout<<std::fixed;
    for(size_t i=0;i<8;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod15(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod15(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<4;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss15(i)<<std::endl;
    }


    std::cout<<"\nGaussKronrod21 \n"<<std::endl;
    for(size_t i=0;i<11;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod21(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod21(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<5;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss21(i)<<std::endl;
    }


    std::cout<<"\nGaussKronrod31 \n"<<std::endl;
    for(size_t i=0;i<16;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod31(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod31(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<8;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss31(i)<<std::endl;
    }


    std::cout<<"\nGaussKronrod41 \n"<<std::endl;
    for(size_t i=0;i<21;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod41(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod41(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<10;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss41(i)<<std::endl;
    }

    std::cout<<"\nGaussKronrod51 \n"<<std::endl;
    for(size_t i=0;i<26;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod51(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod51(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<13;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss51(i)<<std::endl;
    }


    std::cout<<"\nGaussKronrod61 \n"<<std::endl;
    for(size_t i=0;i<31;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::abscissaeGaussKronrod61(i)
            <<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod61(i)<<std::endl;
    }
    std::cout<<std::endl;
    for(size_t i=0;i<15;++i)
    {
        std::cout<<std::setprecision(33)<<QuadratureKronrodValuesType::weightsGauss61(i)<<std::endl;
    }

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=EXIT_SUCCESS;
    test_values();
    return ret;
}
