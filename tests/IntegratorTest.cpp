#include <NIHeaders.h>

#include "gtest/gtest.h"
#include <iostream>
class IntegratorTest : public ::testing::Test
{
protected:
  IntegratorTest() : integrator(200)
  {
  }

  /**
   * This is the closed form of the integral of integrandPeak. The bounds are 0 to 1.
   */
  static double integralPeak(const double alpha);

  /**
   * This is the closed form of the integral of intregrandOscillates. The bounds are 0 to pi.
   */
  static double integralOscillates(const int alpha);

  /**
   * This is the closed form of the integral of integrandInfinite. The bounds are 0 to
   * 40*2^alphaInfinite.
   */
  static float integralInfinite(const float alpha);

  Eigen::Integrator<double> integrator;
};

/**
 * This integrand has a peak of height 4^alphaPeak at x = pi/4.
 */
class IntegrandPeakFunctor
{
public:
  double operator()(const double param) const
  {
    return pow(4., -m_alpha) / (pow(param - M_PI / 4., 2.) + pow(16., -m_alpha));
  }

  /**
   * @param alpha A parameter for varying the peak.
   */
  void setAlpha(const double alpha) {m_alpha = alpha;}

private:
  double m_alpha;
};

/**
 * This integrand oscillates more strongly for increasing alphaOscillates.
 */
class IntegrandOscillatesFunctor
{
public:
  double operator()(const double param) const
  {
    return cos(pow(2., static_cast<double>(m_alpha)) * sin(param));
  }

  /**
   * A parameter for varying the oscillation strength.
   */
  void setAlpha(const double alpha) {m_alpha = alpha;}

private:
  double m_alpha;
};

/**
 * This integrand has an infinite interval. It is well-behaved and tends to 0 quickly.
 */
class IntegrandInfiniteFunctor
{
public:
  float operator()(const float param) const
  {
    return pow(param, 2.) * exp(-param * pow(2, -m_alpha));
  }

  /**
   * A paramater for varying the upper bound.
   */
  void setAlpha(const float alpha) {m_alpha = alpha;}

private:
  float m_alpha;
};

double IntegratorTest::integralPeak(const double alpha)
{
  double factor = pow(4., alpha - 1.);

  return atan((4 - M_PI) * factor) + atan(M_PI * factor);
}

/*
 * The closed form of this integral has the 0th order Bessel function of the first kind as a
 * factor. Literals are used in the absence of a routine for computing this Bessel function.
 */
double IntegratorTest::integralOscillates(const int alpha)
{
  double integrals[11] =
    {
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
    };

  return integrals[alpha];
}

float IntegratorTest::integralInfinite(const float alpha)
{
  float e40 = exp(40.);

  return (e40 - 841.) * pow(2., 3. * alpha + 1.) / e40;
}

template <typename Scalar_>
Scalar_ desiredRelativeError()
{
  return std::numeric_limits<Scalar_>::epsilon() * 50.;
}

template <typename Scalar_>
typename Eigen::Integrator<Scalar_>::QuadratureRule quadratureRules(const size_t i)
{
  static const typename Eigen::Integrator<Scalar_>::QuadratureRule quadratureRules[6] =
    {
      Eigen::Integrator<Scalar_>::GaussKronrod15,
      Eigen::Integrator<Scalar_>::GaussKronrod21,
      Eigen::Integrator<Scalar_>::GaussKronrod31,
      Eigen::Integrator<Scalar_>::GaussKronrod41,
      Eigen::Integrator<Scalar_>::GaussKronrod51,
      Eigen::Integrator<Scalar_>::GaussKronrod61
    };

  return quadratureRules[i];
}

TEST_F(IntegratorTest, qagPeak)
{
  const size_t numKeys = 6;

  IntegrandPeakFunctor integrandPeakFunctor;

  for (size_t i = 0; i < numKeys; ++i)
  {
    Eigen::Integrator<double>::QuadratureRule quadratureRule = quadratureRules<double>(i);

    for (double alpha = 0.; alpha < 18.; ++alpha)
    {
      integrandPeakFunctor.setAlpha(alpha);
      double actual = integrator.quadratureAdaptive(
        integrandPeakFunctor, 0., 1., 0., desiredRelativeError<double>(), quadratureRule);

      double expected = integralPeak(alpha);

      EXPECT_LE(fabs(expected - actual), desiredRelativeError<double>() * fabs(expected));
    }
  }
}

TEST_F(IntegratorTest, qagOscillates)
{
  const size_t numKeys = 6;

  IntegrandOscillatesFunctor integrandOscillatesFunctor;

  for (size_t i = 0; i < numKeys; ++i)
  {
    Eigen::Integrator<double>::QuadratureRule quadratureRule = quadratureRules<double>(i);

    for (int alpha = 0; alpha <= 10; ++alpha)
    {
      integrandOscillatesFunctor.setAlpha(alpha);
      double actual = integrator.quadratureAdaptive(
        integrandOscillatesFunctor, 0., M_PI, 0., desiredRelativeError<double>(), quadratureRule);

      double expected = integralOscillates(alpha);

      EXPECT_LE(fabs(expected - actual), desiredRelativeError<double>() * fabs(expected));
    }
  }
}

TEST_F(IntegratorTest, qagInfinite)
{
  const size_t numKeys = 6;

  IntegrandInfiniteFunctor integrandInfiniteFunctor;

  Eigen::Integrator<float> floatIntegrator(200);

  for (size_t i = 0; i < numKeys; ++i)
  {
    Eigen::Integrator<float>::QuadratureRule quadratureRule = quadratureRules<float>(i);

    for (float alpha = 0.; alpha <= 5.; ++alpha)
    {
      integrandInfiniteFunctor.setAlpha(alpha);
      float actual = floatIntegrator.quadratureAdaptive(
        integrandInfiniteFunctor, 0., 40. * pow(2., alpha), 0., desiredRelativeError<float>(),
        quadratureRule);

      float expected = integralInfinite(alpha);

      EXPECT_LE(fabs(expected - actual), desiredRelativeError<float>() * fabs(expected));
    }
  }
}
