#include <NIHeaders.hpp>

#include <iostream>
#include <iomanip>

// PI must defined for use of quad precision, the GNU C preprocessor value of M_PI
// is double/long double, which will result in reduced accuracy in multiprecision.
#ifndef PI
    #define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899
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


//////////////////////////// Mark's peak function test /////////////////////////
/**
 * This integrand has a peak of height 4^alphaPeak at x = pi/4.
 */
template<typename Scalar>
class IntegrandPeakFunctor
{
public:
  Scalar operator()(const Scalar param) const
  {
    return pow(4., -m_alpha) / (pow(param - PI / 4., 2.) + pow(16., -m_alpha));
  }

  /**
   * @param alpha A parameter for varying the peak.
   */
  void setAlpha(const Scalar alpha) {m_alpha = alpha;}

  static Scalar integralPeak(const Scalar alpha)
  {
    Scalar factor = pow(4., alpha - 1.);
    return atan((4. - PI) * factor) + atan(PI * factor);
  }

private:
  Scalar m_alpha;
};

int test_peak(void)
{
    //typedef float Scalar;
    typedef double Scalar;

    //typedef mpfr::mpreal Scalar;
    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandPeakFunctor<Scalar> IntegrandPeakFunctorType;

    IntegratorType eigenIntegrator(200);
    IntegrandPeakFunctorType integrandPeakFunctor;

    const size_t numKeys = 6;
    for (size_t i = 0; i < numKeys; ++i)
    {
        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(0);

        for (Scalar alpha = 0.; alpha < 18.; ++alpha)
        {
            integrandPeakFunctor.setAlpha(alpha);

            Scalar actual = eigenIntegrator.quadratureAdaptive(integrandPeakFunctor, Scalar(0.),Scalar(1.), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);

            Scalar expected = IntegrandPeakFunctorType::integralPeak(alpha);

            if(fabs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * fabs(expected))
            {
                std::cout << "rule " << i << "\t fabs((Scalar)(expected - actual)) =" << fabs((Scalar)(expected - actual))
                << "\t desiredRelativeError<Scalar>() * fabs(expected)= " << desiredRelativeError<Scalar>() * fabs(expected) << std::endl;
                return EXIT_FAILURE;
            }
        }
    }


    return EXIT_SUCCESS;
}


int main(void)
{
    int ret = EXIT_SUCCESS;
    ret += test_peak();
    return EXIT_SUCCESS;
}
