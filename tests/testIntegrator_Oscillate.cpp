#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

/////////////////////// A generalized example of an oscillatory function ///////////////////
/**
 * This is the closed form of the integral of intregrandOscillates bounded from 0 to pi. 
 * This integrand oscillates more strongly for increasing values of alpha.
 */
template<typename Scalar>
class IntegrandOscillateFunctor
{
public:
    Scalar operator()(const Scalar& param) const
    {
      using std::sin;
      using std::cos;
      using std::pow;
      return cos(pow(2., static_cast<Scalar>(m_alpha)) * sin(param));
    }

    /**
     * \param alpha A parameter for varying the oscillation strength.
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

/**
 * The closed form of this integral has the 0th order Bessel function of the first kind as a
 * factor. Literals are used in the absence of a routine for computing this Bessel function.
 */
template <typename Scalar>
Scalar integralOscillates(const int& alpha)
{
  template <typename Scalar>
  Array<Scalar, 11, 1> integrationValues =
    (Array<Scalar, 11, 1>() <<
       2.4039394306344129,
       0.7033736269566008,
      -1.2476829250428461,
       0.5392569146860977,
      -0.5494616459466271,
       0.4337880026347335,
       0.2908801021737259,
       0.0046251228506773,
      -0.1151503602390470,
      -0.0718346295951386,
       0.0458999248689193
    );

  return integrationValues[alpha];
}

int test_oscillate(void)
{
  std::ofstream fout;
  fout.open("Oscillatory_integration_test_output.txt");

  std::cout<<"Testing "<<std::endl;

  //typedef float Scalar;
  //typedef double Scalar;
  //typedef long double Scalar;
  typedef mpfr::mpreal Scalar;
  Scalar::set_default_prec(16);

  typedef Eigen::Integrator<Scalar> IntegratorType;
  typedef IntegrandOscillateFunctor<Scalar> IntegrandOscillateFunctorType;

  //compute the nodes and weights on the fly
  QuadratureKronrod<Scalar>::computeNodesAndWeights();

  IntegratorType eigenIntegrator(256);
  IntegrandOscillateFunctorType integrandOscillateFunctor;

  bool success = true;
  int counter = 0;
  const size_t numRules = 12;

  for (size_t i = 0; i < numRules; ++i)
  {
    counter = 0;
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

    for (int alpha = 0; alpha < 11; ++alpha)
    {
      success = true;
      integrandOscillateFunctor.setAlpha(alpha);

      Scalar actual = eigenIntegrator.quadratureAdaptive(integrandOscillateFunctor, Scalar(0.), Scalar(M_PI), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
      Scalar expected = integralOscillates(alpha);

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
  ret += test_oscillate();
  return EXIT_SUCCESS;
}