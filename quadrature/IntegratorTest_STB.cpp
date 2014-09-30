#include "Integrator.h"

#include "gtest/gtest.h"
#include <iostream>

//typedef float Scalar; //does not work
typedef double Scalar;
//typedef long double Scalar;

class IntegratorTest : public ::testing::Test
{
public:
    //typedef double Scalar;
protected:
  IntegratorTest() : integrator(200)
  {
  }

  /**
   * This is the closed form of the integral of integrandPeak. The bounds are 0 to 1.
   */
  static Scalar integralPeak(const Scalar alpha)
  {
    Scalar factor = pow(4., alpha - 1.);
    return atan((4 - M_PI) * factor) + atan(M_PI * factor);
}


  Eigen::Integrator<Scalar> integrator;
};

/**
 * This integrand has a peak of height 4^alphaPeak at x = pi/4.
 */
template<typename Scalar>
class IntegrandPeakFunctor
{
public:
  Scalar operator()(const Scalar param) const
  {
    return pow(4., -m_alpha) / (pow(param - M_PI / 4., 2.) + pow(16., -m_alpha));
  }

  /**
   * @param alpha A parameter for varying the peak.
   */
  void setAlpha(const Scalar alpha) {m_alpha = alpha;}

private:
  Scalar m_alpha;
};

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

TEST_F(IntegratorTest, qagPeak)
{

  const size_t numKeys = 6;

  IntegrandPeakFunctor<Scalar> integrandPeakFunctor;

  for (size_t i = 0; i < numKeys; ++i)
  {
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

    for (Scalar alpha = 0.; alpha < 18.; ++alpha)
    {
      integrandPeakFunctor.setAlpha(alpha);
      Scalar actual = integrator.quadratureAdaptive(
        integrandPeakFunctor, 0., 1., 0., desiredRelativeError<Scalar>(), quadratureRule);

      Scalar expected = integralPeak(alpha);

      EXPECT_LE(abs(expected - actual), desiredRelativeError<Scalar>() * abs(expected));
    }
  }
}
