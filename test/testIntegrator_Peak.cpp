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
        // @TODO The usage of constant Pi with fixed precision needs to be changed to the following for multiprecision
        //RealScalar pi = NumTraits<RealScalar>::Pi();
        using std::pow;
        return pow(Scalar(4.), -m_alpha) / (pow(param - Scalar(M_PI) / Scalar(4.), Scalar(2.)) + pow(Scalar(16.), -m_alpha));
    }

    /**
    * \param alpha A parameter for varying the peak.
    */
    void setAlpha(const Scalar& alpha) {m_alpha = alpha;}

    static Scalar integralPeak(const Scalar& alpha)
    {
        using std::pow;
        using std::atan;
        Scalar factor = pow(Scalar(4.), alpha - Scalar(1.));
        return atan((Scalar(4.) - Scalar(M_PI)) * factor) + atan(Scalar(M_PI) * factor);
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

int test_peak(void)
{
    std::ofstream fout;
    fout.open("test/testOutput/Peak_integration_test_output.txt");

    std::cout<<"\nTesting Int [0->1] 4^-alpha/(x-pi/4)^2 + 16^-alpha = atan( (4-pi)4^(alpha-1) )+atan(pi-4^(alpha-1))\n";

    /**
     * When using Multiprecision mpreal types beyond quad precision, it is important to either call
     * computeNodesAndWeights() to calculate nodes and weights on the fly, or to repopulate the
     * QuadratureKronrod.h tabulated array values by first changing the setprecision and output 
     * precision values in outputKronrodNodesWeights.cpp to greater than the value of precision
     * required for the subseuent integration calculations desired.
     */
     
    //typedef float Scalar;
    typedef double Scalar;
    //typedef long double Scalar;
    
    /**
     * typedef mpfr::mpreal Scalar;
     * Scalar::set_default_prec(6);
     * QuadratureKronrod<Scalar>::computeNodesAndWeights();
     */

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandPeakFunctor<Scalar> IntegrandPeakFunctorType;

    IntegratorType eigenIntegrator(100);
    IntegrandPeakFunctorType integrandPeakFunctor;

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
            integrandPeakFunctor.setAlpha(alpha);

            Scalar actual = eigenIntegrator.quadratureAdaptive(integrandPeakFunctor, Scalar(0.),Scalar(1.), Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);
            Scalar expected = IntegrandPeakFunctorType::integralPeak(alpha);

            using std::abs;
            if(abs((Scalar)(expected - actual)) > desiredRelativeError<Scalar>() * abs(expected) 
                || isnan(abs((Scalar)(expected - actual))))
            {
                fout << "\nrule " << i << "\n abs(expected - actual) = " << abs(expected - actual)
                     << "\n desiredRelativeError<Scalar>() * abs(expected) = "
                     << desiredRelativeError<Scalar>() * abs(expected)<<std::endl;

                fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
                success = false;
            }
            else
            {
                fout << "\nrule " << i << "\n abs(expected - actual) = " << abs(expected - actual)
                     << "\n desiredRelativeError<Scalar>() * Abs(expected) = "
                     << desiredRelativeError<Scalar>() * abs(expected) << std::endl;
                          
                fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
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
    ret += test_peak();
    return EXIT_SUCCESS;
}
