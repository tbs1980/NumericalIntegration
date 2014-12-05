#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

////////////////////////// example from Quadpackcpp ///////////////////////////
/**
 * Int [0->1] x^alog(1/x) = 1/(a+1)^2
 */
template<typename Scalar>
class IntegrandLogPowFunctor
{
public:
  Scalar operator()(const Scalar& param) const
  {
    using std::pow;
    using std::log;
    return pow(param, m_alpha) * log(1/param);
  }

  /**
   * \param alpha A parameter for varying the upper bound.
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

int test_logpow(void)
{
  std::ofstream fout;
  fout.open("LogPow_integration_test_output.txt");

  std::cout<<"\nTesting Int [0->1] x^a*log(1/x) = 1/(a+1)^2\n";

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
  IntegrandLogPowFunctorType integrandLogPowFunctor;

  bool success = true;
  int counter = 0;
  const size_t numRules = 12;

  for (size_t i = 0; i < numRules; ++i)
  {
  	counter = 0;
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

    for (Scalar alpha = 0.; alpha < 18.; ++alpha)
    {
      success = true;
      integrandLogPowFunctor.setAlpha(alpha);

      Scalar actual = eigenIntegrator.quadratureAdaptive(integrandLogPowFunctor, Scalar(0.),Scalar(1.), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
      Scalar expected = IntegrandLogPowFunctorType::exact_value_in_01(alpha);

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
  ret += test_logpow();
  return EXIT_SUCCESS;
}
