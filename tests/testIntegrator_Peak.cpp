#include <NIHeaders.hpp>

#include <iostream>
#include <iomanip>

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

            if(Abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * Abs(expected))
            {
                std::cout<<"rule "<<i<<"\t abs((Scalar)(expected - actual)) ="<<Abs((Scalar)(expected - actual))
                <<"\t desiredRelativeError<Scalar>() * Abs(expected)= "<<desiredRelativeError<Scalar>() * Abs(expected)<<std::endl;
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
