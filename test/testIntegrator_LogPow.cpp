/**
 * \file testIntegrator_LogPow.cpp
 * This file is a unit test for Integrator.h and its' associated files.
 * The test function is a decaying function which begins at positive 
 * infinity and appraoches the final asymptotic value rapidly.
 */

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
        using std::log;
        using std::pow;

        return pow(param, m_alpha) * log(1/param);
    }

    /**
    * \param alpha A parameter for varying the upper bound.
    */
    void setAlpha(const Scalar& alpha)
    {
        m_alpha = alpha;
    }

    static Scalar exact_value_in_01(const Scalar& alpha)
    {
        Scalar a1 = alpha + 1;
        return 1/(a1*a1);
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

int test_logpow(void)
{
    using std::abs;
    using std::isnan;
    
    std::ofstream fout;
    fout.open("test/testOutput/LogPow_integration_test_output.txt");

    std::cout<<"\nTesting Int [0->1] x^a*log(1/x) = 1/(a+1)^2\n";
     
    // typedef float Scalar;
    typedef double Scalar;
    // typedef long double Scalar;
    // typedef mpfr::mpreal Scalar;
    // Scalar::set_default_prec(500);   // \detail This sets the number of bits of precision; each signficant figure desired will require 4 bits.
    // QuadratureKronrod<Scalar>::computeNodesAndWeights(); // \detail Utilizing precision beyond double requires nodes to be computed at runtime, because of the manner that the static values are truncated when they are assigned at compile time.

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandLogPowFunctor<Scalar> IntegrandLogPowFunctorType;

    IntegratorType eigenIntegrator(10000);  // \detail The number of subintervals must be increased by more than 100X the precision requested.
    IntegrandLogPowFunctorType integrandLogPowFunctor;

    bool success = true;
    const Scalar alphaLimit = Scalar(18.);
    const Index numRules = 12;

    for (Scalar alpha = Scalar(0.); alpha < alphaLimit; ++alpha)
    {
        success = true;
        integrandLogPowFunctor.setAlpha(alpha);
        
        for (Index i = 0; i < numRules; ++i)
        {
            Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);
            
            Scalar actual = eigenIntegrator.quadratureAdaptive(integrandLogPowFunctor, Scalar(0.),Scalar(1.), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
            Scalar expected = IntegrandLogPowFunctorType::exact_value_in_01(alpha);

            if(abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * abs(expected) 
                || isnan(abs((Scalar)(expected - actual))))
            {
                success = false;

                if(i == numRules-1)
                {
                    fout << "\nPeak Test could not pass Alpha = " << alpha
                         << "\nrule " << i+1 << "\n abs(expected - actual) = " << abs(expected - actual)
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
                fout << "\nrule " << i+1 << "\n abs(expected - actual) = " << abs(expected - actual)
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

    if(success)
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
    ret += test_logpow();
    return EXIT_SUCCESS;
}
