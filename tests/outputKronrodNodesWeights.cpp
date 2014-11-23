#include <NIHeaders.h>

#include <iostream>
#include <fstream>
#include <iomanip>

int test_values()
{
    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);
    typedef Eigen::QuadratureKronrod<Scalar> QuadratureKronrodValuesType;

    QuadratureKronrodValuesType::ComputeNodesAndWeights();

    int outputDigits = 100;
    std::cout<<std::fixed;
/*
    //--------15--------//
    std::cout<<"\nGaussKronrod15 \n"<<std::endl;
    for(size_t i=0;i<8;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod15(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<8;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod15(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<4;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss15(i)<<","<<std::endl;
    }

    //--------21--------//
    std::cout<<"\nGaussKronrod21 \n"<<std::endl;
    for(size_t i=0;i<11;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod21(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<11;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod21(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<5;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss21(i)<<","<<std::endl;
    }

    //--------31--------//
    std::cout<<"\nGaussKronrod31 \n"<<std::endl;
    for(size_t i=0;i<16;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod31(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<16;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod31(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<8;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss31(i)<<","<<std::endl;
    }

    //--------41--------//
    std::cout<<"\nGaussKronrod41 \n"<<std::endl;
    for(size_t i=0;i<21;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod41(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<21;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod41(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<10;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss41(i)<<","<<std::endl;
    }

    //--------51--------//
    std::cout<<"\nGaussKronrod51 \n"<<std::endl;
    for(size_t i=0;i<26;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod51(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<26;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod51(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<13;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss51(i)<<","<<std::endl;
    }

    //--------61--------//
    std::cout<<"\nGaussKronrod61 \n"<<std::endl;
    for(size_t i=0;i<31;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod61(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<31;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod61(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<15;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss61(i)<<","<<std::endl;
    }

*/
    //--------71--------//
    std::cout<<"\nGaussKronrod71 \n"<<std::endl;
    for(size_t i=0;i<36;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod71(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<36;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod71(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<18;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss71(i)<<","<<std::endl;
    }

    //--------81--------//
    std::cout<<"\nGaussKronrod81 \n"<<std::endl;
    for(size_t i=0;i<41;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod81(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<41;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod81(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<20;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss81(i)<<","<<std::endl;
    }

    //--------91--------//
    std::cout<<"\nGaussKronrod91 \n"<<std::endl;
    for(size_t i=0;i<46;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod91(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<46;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod91(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<23;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss91(i)<<","<<std::endl;
    }
/*
    //--------101--------//
    std::cout<<"\nGaussKronrod101 \n"<<std::endl;
    for(size_t i=0;i<51;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod101(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<51;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod101(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<25;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss101(i)<<","<<std::endl;
    }

    //--------201--------//
    std::cout<<"\nGaussKronrod201 \n"<<std::endl;
    for(size_t i=0;i<101;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod201(i)<<","<<std::endl;
    }
    std::cout<<std::endl;
    
    for(size_t i=0;i<101;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod201(i)<<","<<std::endl;
    }
    std::cout<<std::endl;

    for(size_t i=0;i<50;++i)
    {
        std::cout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss201(i)<<","<<std::endl;
    }
*/
    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=EXIT_SUCCESS;
    test_values();
    return ret;
}
