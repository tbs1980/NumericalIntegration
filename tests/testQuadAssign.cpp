#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>


template <typename RealType>
class QuadratureKronrod
{
public:
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef typename LaurieGautschiPolicy::IndexType IndexType;
    QuadratureKronrod()
    {
        LaurieGautschiPolicy::mpkonrad15(abscissaeGaussKronrod15,weightsGaussKronrod15);
        LaurieGautschiPolicy::mpkonrad21(abscissaeGaussKronrod21,weightsGaussKronrod21);
        LaurieGautschiPolicy::mpkonrad31(abscissaeGaussKronrod31,weightsGaussKronrod31);
        LaurieGautschiPolicy::mpkonrad41(abscissaeGaussKronrod41,weightsGaussKronrod41);
        LaurieGautschiPolicy::mpkonrad51(abscissaeGaussKronrod51,weightsGaussKronrod51);
        LaurieGautschiPolicy::mpkonrad61(abscissaeGaussKronrod61,weightsGaussKronrod61);
    }
  Array<RealType, 8, 1> abscissaeGaussKronrod15;
  Array<RealType, 8, 1> weightsGaussKronrod15;
  Array<RealType, 4, 1> weightsGauss15;

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

};



int main(void)
{
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);

    typedef QuadratureKronrod<RealType> QuadratureKronrodType;
    typedef QuadratureKronrodType::IndexType IndexType;

    QuadratureKronrodType qr;

    std::cout<<"\n GaussKronrod15\n"<<std::endl;
    std::cout<<std::fixed;
    for(IndexType i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<qr.abscissaeGaussKronrod15(i)<<"\t"<<qr.weightsGaussKronrod15(i)<<std::endl;
    }

    return EXIT_SUCCESS;
}
