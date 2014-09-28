#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>


template <typename RealType>
class QuadratureKronrod
{
public:

    //typedef typename LaurieGautschiPolicy::IndexType IndexType;
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
  static const Array<RealType, 8, 1> abscissaeGaussKronrod15;
  static const Array<RealType, 8, 1> weightsGaussKronrod15;
  static const Array<RealType, 4, 1> weightsGauss15;

    /*
  Array<RealType, 11, 1> abscissaeGaussKronrod21;
  Array<RealType, 11, 1> weightsGaussKronrod21;
  Array<RealType, 5, 1> weightsGauss21;

  Array<RealType, 16, 1> abscissaeGaussKronrod31;
  Array<RealType, 16, 1> weightsGaussKronrod31;
  Array<RealType, 8, 1> weightsGauss31;

  Array<RealType, 21, 1> abscissaeGaussKronrod41;
  Array<RealType, 21, 1> weightsGaussKronrod41;
  Array<RealType, 10, 1> weightsGauss41;

  Array<RealType, 26, 1> abscissaeGaussKronrod51;
  Array<RealType, 26, 1> weightsGaussKronrod51;
  Array<RealType, 13, 1> weightsGauss51;

  Array<RealType, 31, 1> abscissaeGaussKronrod61;
  Array<RealType, 31, 1> weightsGaussKronrod61;
  Array<RealType, 15, 1> weightsGauss61;
  */
};

template <typename RealType>
const Array<RealType, 8, 1> QuadratureKronrod<RealType>::abscissaeGaussKronrod15 =
    QuadratureKronrod<RealType>::LaurieGautschiPolicy::mpkonrad15abscissae();

template <typename RealType>
const Array<RealType, 8, 1> QuadratureKronrod<RealType>::weightsGaussKronrod15 =
    QuadratureKronrod<RealType>::LaurieGautschiPolicy::mpkonrad15weights();

int main(void)
{
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    //typedef double RealType;

    typedef QuadratureKronrod<RealType> QuadratureKronrodType;
    typedef QuadratureKronrodType::LaurieGautschiPolicy LaurieGautschiPolicy;
    //typedef QuadratureKronrodType::IndexType IndexType;


    std::cout<<"\n GaussKronrod15\n"<<std::endl;
    std::cout<<std::fixed;
    for(int i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<QuadratureKronrodType::abscissaeGaussKronrod15(i)<<"\t"
            <<QuadratureKronrodType::weightsGaussKronrod15(i)<<std::endl;
    }

    /*
    std::cout<<"\n try 2 GaussKronrod15\n"<<std::endl;
    Array<RealType, 8, 1> xout;
    Array<RealType, 8, 1> wout;
    xout = LaurieGautschiPolicy::mpkonrad15abscissae();
    wout = LaurieGautschiPolicy::mpkonrad15weights();

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
