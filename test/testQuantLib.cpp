#include <NumericalIntegration.h>
#include <ql/math/integrals/kronrodintegral.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>

/**
 * \param epsilon Relative machine precision.
 */
template <typename Scalar>
Scalar desiredRelativeError()
{
  return Eigen::NumTraits<Scalar>::epsilon() * 50.;
}


template<typename T>
T mySin(T x)
{
    return std::sin(x);
}

template<typename Scalar>
class IntegrandSineFunctor
{
public:
    Scalar operator()(const Scalar& param) const
    {
      using std::sin;
      return sin(param);
    }
};


boost::function<double (double x)> f;

void QuantLibSineIntegration(void)
{
    typedef QuantLib::Real Scalar;
    f = &mySin<Scalar>;
    QuantLib::Size maxEvaluations = 2000;
    Scalar tolerance = desiredRelativeError<QuantLib::Real>();//1.0e-6;

    QuantLib::GaussKronrodAdaptive I(tolerance,maxEvaluations);

    Scalar a(0);
    Scalar b(M_PI);

    Scalar expected(2.);
    Scalar calculated = I(f,a,b);

    std::cout<<"|calculated - expected| from QuantLib"<<std::endl;
    std::cout<<std::setprecision(10)<<std::abs(calculated - expected)<<"\n"<<std::endl;
}


void EigenSineIntegration()
{
    typedef QuantLib::Real Scalar;
    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandSineFunctor<Scalar> IntegrandSineFunctorType;

    IntegratorType eigenIntegrator(1000);
    IntegrandSineFunctorType integrandSineFunctor;

    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod15;

    Scalar expected(2);
    Scalar calculated = eigenIntegrator.quadratureAdaptive(integrandSineFunctor,
        Scalar(0.), Scalar(M_PI), Scalar(0.), desiredRelativeError<Scalar>(),
        quadratureRule);

    std::cout<<"|calculated - expected| from Eigen"<<std::endl;
    std::cout<<std::scientific<<std::setprecision(10)<<std::abs(calculated - expected)<<"\n"<<std::endl;
}

int main(void)
{
    QuantLibSineIntegration();
    EigenSineIntegration();
    return 0;
}
