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
    std::cout<<"Testing Int [0->1] 4^-alpah/(x-pi/4)^2 + 16^-alpha = atan( (4-pi)4^(alpha-1) )+atan(pi-4^(alpha-1))"<<std::endl;
    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(53);

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandPeakFunctor<Scalar> IntegrandPeakFunctorType;

    //compute the nodes and weights on the fly
    QuadratureKronrod<Scalar>::ComputeNodesAndWeights();

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
                std::cout << "\nrule " << i << "\n fabs(expected - actual) =" << fabs(expected - actual)
                          << "\n desiredRelativeError<Scalar>() * Abs(expected)= "
                          << desiredRelativeError<Scalar>() * fabs(expected)<<std::endl;

                std::cout << "erroCode =" << eigenIntegrator.errorCode() << std::endl;

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
