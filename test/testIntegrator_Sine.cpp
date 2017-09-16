/**
 * \file testIntegrator_Sine.cpp
 * This file is a unit test for Integrator.h and its' associated files.
 * The test is an (oscillatory) Sine function from zero to pi.
 */

#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;

///////////////////////////// A simple example of sin(x) ///////////////////////
template<typename Scalar>
class IntegrandSineFunctor
{
public:
    Scalar operator()(const Scalar& param) const
    {
        using std::sin;
        return sin(param);
    }
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

int test_sine(void)
{
    using std::abs;
    using std::isnan;
    
    std::ofstream fout;
    fout.open("test/testOutput/Sine_integration_test_output.txt");

    std::cout<<"\nTesting Int [0->Pi] sin(x) = 2\n";

    // typedef float Scalar;
    typedef double Scalar;
    // typedef long double Scalar;
    // typedef mpfr::mpreal Scalar;    // \detail Performing this test using multiprecision requires changing from M-PI to NumTraits<Scalar>::PI();
    // Scalar::set_default_prec(500);  // \detail This sets the number of bits of precision; each signficant figure desired will require 4 bits.
    // QuadratureKronrod<Scalar>::computeNodesAndWeights(); // \detail Utilizing precision beyond long double requires nodes to be computed at runtime, because of the manner that the static values are truncated when they are assigned at compile time.

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandSineFunctor<Scalar> IntegrandSineFunctorType;

    IntegratorType eigenIntegrator(1000); // \detail The number of subintervals must be increased by more than 100X the precision requested.
    IntegrandSineFunctorType integrandSineFunctor;

    bool success = true;
    const Index numRules = 12;

    for (Index i = 0; i < numRules; ++i)
    {
        success = true;

        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

        Scalar actual = eigenIntegrator.quadratureAdaptive(integrandSineFunctor, Scalar(0.), Scalar(M_PI), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
        // \detail The usage of NumTraits<Scalar>::Pi() is required for multiprecision
        // Scalar actual = eigenIntegrator.quadratureAdaptive(integrandSineFunctor, Scalar(0.), NumTraits<Scalar>::Pi(), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
        Scalar expected = Scalar(2.);

        if (abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * abs(expected) ||
            isnan((Scalar)(expected - actual)))
        {
            fout << "\nrule " << i+1 << "\n abs(expected - actual) = " << abs(expected - actual)
                 << "\n desiredRelativeError<Scalar>() * abs(expected) = "
                 << desiredRelativeError<Scalar>() * abs(expected) << std::endl;

          fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
          success = false;
        }
        else
        {
            fout << "\nrule " << i+1 << "\n abs(expected - actual) = " << abs(expected - actual)
                 << "\n desiredRelativeError<Scalar>() * abs(expected) = "
                 << desiredRelativeError<Scalar>() * abs(expected) << std::endl;

            fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
            fout << "Success!\n ";
        }

        if(success)
        {
            fout << "\n\tTest Succeeded!\n" << std::endl;
            fout.close();
            break;
        }
        else
        {
            fout <<"\n\tTest Failed.\n" << std::endl;
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
        std::cout << std::endl << "\tTest Failed.\n" << std::endl;
        return EXIT_FAILURE;
    }
}

int main(void)
{
    int ret = EXIT_SUCCESS;
    ret += test_sine();
    return EXIT_SUCCESS;
}
