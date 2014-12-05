#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

////////////////////////// example from Quadpackcpp ///////////////////////////
/**
 * This integrand has an infinite interval. It is well-behaved and tends to 0 quickly.
 */
template<typename Scalar>
class IntegrandPowFunctor
{
public:
  Scalar operator()(const Scalar& param) const
  {
    using std::pow;
    using std::exp;
    return pow(param, 2.) * exp(-param * pow(2, -m_alpha));
  }

  /**
   * A paramater for varying the upper bound.
   */
  void setAlpha(const Scalar& alpha) {m_alpha = alpha;}

  static Scalar exact_value_in_01(const Scalar& alpha)
  {
    Scalar a1 = alpha+1.;
    return 1./(a1*a1);
  }

private:
  Scalar m_alpha;
};

/**
 * \param epsilon Relative machine precision.
 */
template <typename Scalar>
Scalar desiredRelativeError()
{
  return Eigen::NumTraits<Scalar>::epsilon() * 50.;
}

template <typename Scalar>
typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules(const size_t& i)
{
  static const typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules[12] =
    {
      Eigen::Integrator<Scalar>::GaussKronrod15,
      Eigen::Integrator<Scalar>::GaussKronrod21,
      Eigen::Integrator<Scalar>::GaussKronrod31,
      Eigen::Integrator<Scalar>::GaussKronrod41,
      Eigen::Integrator<Scalar>::GaussKronrod51,
      Eigen::Integrator<Scalar>::GaussKronrod61,
      Eigen::Integrator<Scalar>::GaussKronrod71,
      Eigen::Integrator<Scalar>::GaussKronrod81,
      Eigen::Integrator<Scalar>::GaussKronrod91,
      Eigen::Integrator<Scalar>::GaussKronrod101,
      Eigen::Integrator<Scalar>::GaussKronrod121,
      Eigen::Integrator<Scalar>::GaussKronrod201
    };

  return quadratureRules[i];
}

Scalar IntegratorTest::integralPow(const Scalar& alpha)
{
  Scalar e40 = exp(40.);
  using std::pow;
  return (e40 - 841.) * pow(2., 3. * alpha + 1.) / e40;
}

int test_pow(void)
{
    std::ofstream fout;
    fout.open("Pow_integration_test_output.txt");

    std::cout<<"Testing "<<std::endl;

    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(114);

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandLogPowFunctor<Scalar> IntegrandLogPowFunctorType;

    //compute the nodes and weights on the fly
    QuadratureKronrod<Scalar>::computeNodesAndWeights();

    IntegratorType eigenIntegrator(256);
    IntegrandInfiniteFunctorType integrandInfiniteFunctor;

    bool success = true;
    int counter = 0;
    const size_t numKeys = 12;

    for (size_t i = 0; i < numKeys; ++i)
    {
      counter = 0;
        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

        for (Scalar alpha = 0.; alpha < 18.; ++alpha)
        {
            success = true;
            integrandInfiniteFunctor.setAlpha(alpha);

            using std::pow;
            Scalar actual = floatIntegrator.quadratureAdaptive(integrandPowFunctor, Scalar(0.), Scalar(40. * pow(2., alpha)), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);

            Scalar expected = IntegrandPowFunctorType::exact_value_in_01(alpha);

            using std::abs;
            if(abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * abs(expected))
            {
                fout << "\nrule " << i << "\n Abs(expected - actual) =" << abs(expected - actual)
                          << "\n desiredRelativeError<Scalar>() * Abs(expected)= "
                          << desiredRelativeError<Scalar>() * abs(expected) << std::endl;

                fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
                success = false;
            }
            else
            {
                fout << "\nrule " << i << "\n abs(expected - actual) =" << abs(expected - actual)
                          << "\n desiredRelativeError<Scalar>() * abs(expected)= "
                          << desiredRelativeError<Scalar>() * abs(expected) << std::endl;
                          
                fout << "Success!\n";
                counter++;
            }
        }

        if(success && counter == 18)    
        {
          fout << "\n  Test Succeeded!\n" << std::endl;
          fout.close();
          break;
        }
        else
        {
          fout <<"\n  Test Failed.\n" << std::endl;
        }
    }

    fout.close();

    if(success && counter == 18)
    {
      std::cout << std::endl << "  Test Succeeded!\n" << std::endl;
      return EXIT_SUCCESS;
    }
    else
    {
      std::cout << std::endl << "  Test Failed.\n" << std::endl;
      return EXIT_FAILURE;
    }
}

int main(void)
{
    int ret = EXIT_SUCCESS;
    ret += test_pow();
    return EXIT_SUCCESS;
}