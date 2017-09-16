/**
 * \file testIntegrator_Peak.cpp
 * This file is a unit test for Integrator.h and its' associated files.
 * The test function peak varies with the parameter alpha to create an
 * increasingly sharp point at the abscissa equal to pi/4.
 */

 #include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

/**
 * This integrand has a peak of height 4^alphaPeak at x = pi/4.
 */
template<typename Scalar>
class IntegrandPeakFunctor
{
public:
    Scalar operator()(const Scalar& param) const
    {
        using std::pow;

        return pow(Scalar(4.), -m_alpha) / (pow(param-Scalar(M_PI)/Scalar(4.), Scalar(2.)) + pow(Scalar(16.), -m_alpha));
        // \detail The usage of NumTraits<Scalar>::Pi() is required for multiprecision
        // return pow(Scalar(4.), -m_alpha) / (pow(param-NumTraits<Scalar>::Pi() / Scalar(4.), Scalar(2.)) + pow(Scalar(16.), -m_alpha));
    }

    /**
    * \param alpha A parameter for varying the peak.
    */
    void setAlpha(const Scalar& alpha)
    {
        m_alpha = alpha;
    }

    static Scalar integralPeak(const Scalar& alpha)
    {
        using std::atan;
        using std::pow;

        return atan((Scalar(4.) - Scalar(M_PI))*pow(Scalar(4.), alpha - Scalar(1.))) + atan(Scalar(M_PI)*pow(Scalar(4.), alpha - Scalar(1.)));
        // \detail The usage of NumTraits<Scalar>::Pi() is required for multiprecision
        // return atan((Scalar(4.) - NumTraits<Scalar>::Pi())*pow(Scalar(4.), alpha - Scalar(1.))) + atan(NumTraits<Scalar>::Pi()*pow(Scalar(4.), alpha - Scalar(1.)));
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
    return NumTraits<Scalar>::epsilon() * Scalar(50.);
}

template <typename Scalar>
typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules(const Index& i)
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

int test_peak(void)
{
    using std::abs;
    using std::isnan;
    
    std::ofstream fout;
    fout.open("test/testOutput/Peak_integration_test_output.txt");

    std::cout<<"\nTesting Int [0->1] 4^-alpha/((x-pi/4)^2 + 16^-alpha) = atan((4-pi)*4^(alpha-1)) + atan(pi*4^(alpha-1))\n";

    // typedef float Scalar;            // \details float precision will not pass beyond alphaLimit = 8.
    typedef double Scalar;              // \details double precision will not pass beyond alphaLimit = 8.
    // typedef long double Scalar;      // \details long double precision will not pass beyond alphaLimit = 10.
    // typedef mpfr::mpreal Scalar;     // \detail Performing this test using multiprecision requires changing from M_PI to NumTraits<Scalar>::PI();
    // Scalar::set_default_prec(350);   // \detail This sets the number of bits of precision; each signficant figure desired will require 4 bits.
    // QuadratureKronrod<Scalar>::computeNodesAndWeights(); // \detail Utilizing precision beyond long double requires nodes to be computed at runtime, because of the manner that the static values are truncated when they are assigned at compile time.
    

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandPeakFunctor<Scalar> IntegrandPeakFunctorType;

    IntegratorType eigenIntegrator(100000);  // \detail The number of subintervals must be increased to roughly 100X the precision requested.
    IntegrandPeakFunctorType integrandPeakFunctor;

    bool success = true;
    const Scalar alphaLimit = Scalar(15.);
    const Index numRules = 12;

    for (Scalar alpha = Scalar(0.); alpha < alphaLimit; ++alpha)
    {
        success = true;
        integrandPeakFunctor.setAlpha(alpha);

        for (Index i = 0; i < numRules; ++i)
        {
            Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

            Scalar actual = eigenIntegrator.quadratureAdaptive(integrandPeakFunctor, Scalar(0.), Scalar(1.), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
            Scalar expected = IntegrandPeakFunctorType::integralPeak(alpha);

            if (abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * abs((Scalar)expected) 
                || isnan(abs((Scalar)(expected - actual))))
            {
                success = false;

                if (i == numRules-1)
                {
                    fout << "\nPeak Test could not pass Alpha = " << alpha
                         << "\nrule " << i << "\n abs(expected - actual) = " << abs(expected - actual)
                         << "\n desiredRelativeError<Scalar>() * abs(expected) = "
                         << desiredRelativeError<Scalar>() * abs(expected) << std::endl;
                          
                    fout << "errorCode = " << eigenIntegrator.errorCode() << "\nTest aborted after Fail";

                    std::cout << "\nTest aborted after failing for Alpha = " << alpha
                              << "\n\tTest Failed.\n" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            else
            {
                fout << "\nrule " << i << "\n abs(expected - actual) = " << abs(expected - actual)
                     << "\n desiredRelativeError<Scalar>() * abs(expected) = "
                     << desiredRelativeError<Scalar>() * abs(expected) << std::endl;
                          
                fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
                fout << "alpha = " << alpha << std::endl;
                fout << "Success!\n";
                success = true;
                break;
            }
        }
    }

    fout.close();

    if (success)
    {
        std::cout << std::endl << "\tTest Succeeded!\n" << std::endl;
        return EXIT_SUCCESS;
    }
    else
    {
        std::cout << "\tTest Failed.\n" << std::endl;
        return EXIT_FAILURE;
    }
}


int main(void)
{
    int ret = EXIT_SUCCESS;
    ret += test_peak();
    return EXIT_SUCCESS;
}
