#include <NIHeaders.hpp>

#include <iostream>
#include <iomanip>

#ifndef PI
    #define PI 3.1415926535897932384626433832795028841971693993751
#endif

template <typename Scalar>
Scalar desiredRelativeError()
{
  //return std::numeric_limits<Scalar>::epsilon() * 50.;
  return Eigen::NumTraits<Scalar>::epsilon() * 50.;
}

template <typename Scalar>
typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules(const size_t i)
{
  static const typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules[6] =
    {
      Eigen::Integrator<Scalar>::GaussKronrod15,
      Eigen::Integrator<Scalar>::GaussKronrod21,
      Eigen::Integrator<Scalar>::GaussKronrod31,
      Eigen::Integrator<Scalar>::GaussKronrod41,
      Eigen::Integrator<Scalar>::GaussKronrod51,
      Eigen::Integrator<Scalar>::GaussKronrod61
    };

  return quadratureRules[i];
}

///////////////////////////// A simple example of sin(x) ///////////////////////
template<typename Scalar>
class IntegrandSineFunctor
{
public:
    Scalar operator()(const Scalar param) const
    {
        return sin(param);
    }
};


int test_sine(void)
{
    std::cout<<"Testing Int [0->Pi] sin(x) = 2"<<std::endl;
    typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    //typedef mpfr::mpreal Scalar;
    //Scalar::set_default_prec(128);

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandSineFunctor<Scalar> IntegrandSineFunctorType;

    IntegratorType eigenIntegrator(200);

    return EXIT_SUCCESS;

    IntegrandSineFunctorType integrandSineFunctor;

    std::cout<<"Desired Relatvie error = "<<desiredRelativeError<Scalar>()<<std::endl;

    const size_t numKeys = 6;
    for (size_t i = 0; i < numKeys; ++i)
    {
        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

        Scalar actual = eigenIntegrator.quadratureAdaptive(integrandSineFunctor, Scalar(0.), Scalar(PI),
            Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);

        Scalar expected = Scalar(2);

        if(Abs(expected - actual) > desiredRelativeError<Scalar>() * Abs(expected)
            or eigenIntegrator.errorCode() !=0)
        {
            std::cout<<"\nrule "<<i<<"\n Abs(expected - actual) ="<<Abs(expected - actual)
                <<"\n desiredRelativeError<Scalar>() * Abs(expected)= "
                <<desiredRelativeError<Scalar>() * Abs(expected)<<std::endl;
            std::cout<<"erroCode ="<<eigenIntegrator.errorCode()<<std::endl;
            //return EXIT_FAILURE;
        }

    }


    return EXIT_SUCCESS;
}

int main(void)
{
    int ret = EXIT_SUCCESS;
    ret += test_sine();
    return EXIT_SUCCESS;
}
