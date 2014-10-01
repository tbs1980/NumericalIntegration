#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

template <typename RealType>
class QuadratureKronrodTest
{
public:

    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;

    QuadratureKronrodTest()
    {
        if(s_ready==false)
        {
            std::cout<<"computing"<<std::endl;
            s_abscissaeGaussKronrod15 = LaurieGautschiPolicy::mpkonrad15abscissae();
            s_weightsGaussKronrod15 = LaurieGautschiPolicy::mpkonrad15weights();
            s_ready = true;
        }
        else
        {
            std::cout<<"already computed"<<std::endl;
        }

    }

    Array<RealType, 8, 1> abscissaeGaussKronrod15()
    {
        return s_abscissaeGaussKronrod15;
    }

    Array<RealType, 8, 1> weightsGaussKronrod15()
    {
        return s_weightsGaussKronrod15;
    }

//private:

    static Array<RealType, 8, 1> s_abscissaeGaussKronrod15 ;
    static Array<RealType, 8, 1> s_weightsGaussKronrod15;
    static Array<RealType, 4, 1> s_weightsGauss15;
    static bool s_ready;

};

template <typename RealType>
bool QuadratureKronrodTest<RealType>::s_ready=false;

template <typename RealType>
Array<RealType, 8, 1> QuadratureKronrodTest<RealType>::s_abscissaeGaussKronrod15 =
    QuadratureKronrodTest<RealType>::LaurieGautschiPolicy::mpkonrad15abscissae();

template <typename RealType>
Array<RealType, 8, 1> QuadratureKronrodTest<RealType>::s_weightsGaussKronrod15 =
    QuadratureKronrodTest<RealType>::LaurieGautschiPolicy::mpkonrad15weights();


int main(void)
{
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);



    typedef QuadratureKronrodTest<RealType> QuadratureKronrodType;
    typedef QuadratureKronrodType::LaurieGautschiPolicy LaurieGautschiPolicy;


    std::cout<<"\n GaussKronrod15-0\n"<<std::endl;
    std::cout<<std::fixed;
    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<QuadratureKronrodType::s_abscissaeGaussKronrod15(i)<<"\t"
            <<QuadratureKronrodType::s_weightsGaussKronrod15(i)<<std::endl;
    }

    QuadratureKronrodType qk;


    std::cout<<"\n GaussKronrod15\n"<<std::endl;
    std::cout<<std::fixed;
    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<qk.abscissaeGaussKronrod15()(i)<<"\t"
            <<qk.weightsGaussKronrod15()(i)<<std::endl;
    }

    QuadratureKronrodType qk1;


    std::cout<<"\n GaussKronrod15-1\n"<<std::endl;
    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<qk1.abscissaeGaussKronrod15()(i)<<"\t"
            <<qk1.weightsGaussKronrod15()(i)<<std::endl;
    }

    /*
    std::cout<<"\n try 2 GaussKronrod15\n"<<std::endl;
    Array<RealType, 8, 1> xout;
    Array<RealType, 8, 1> wout;
    xout = LaurieGautschiPolicy::mpkonrad15abscissae();
    wout = LaurieGautschiPolicy::mpkonrad15weights();

    std::cout<<std::fixed;
    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<xout(i)<<"\t"
            <<wout(i)<<std::endl;
    }

    std::cout<<"\n try 3 GaussKronrod15\n"<<std::endl;
    Array<RealType, 8, 1> xout1;
    Array<RealType, 8, 1> wout1;
    //xout = LaurieGautschiPolicy::mpkonrad15abscissae();
    //wout = LaurieGautschiPolicy::mpkonrad15weights();
    LaurieGautschiPolicy::mpkonrad15(xout1,wout1);

    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<xout1(i)<<"\t"
            <<wout1(i)<<std::endl;
    }
    */

    return EXIT_SUCCESS;
}
