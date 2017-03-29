/**
 *\file Integration.cpp
 * An example illustrating the use of numerical integration module in Eigen.
 */

#include <NumericalIntegration.h>

#include <iostream>
#include <iomanip>

/** 
 *  We consider the example from:
 *
 *       http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html
 *
 *       int_0^1 x^{-1/2} log(x) dx = -4
 *
 *  The integrator expects the user to provide a functor as shown below.
 */

template<typename Scalar>
class IntegrandExampleFunctor
{
public:
    IntegrandExampleFunctor(const Scalar alpha)
        : m_alpha(alpha)
    {
        assert(alpha>0);
    }

    Scalar operator()(const Scalar x) const
    {
        assert(x>0);
        return log(m_alpha*x) / sqrt(x);
    }

    void setAlpha(const Scalar alpha)
    {
        m_alpha = alpha;
    }

private:
    Scalar m_alpha;
};

int main(void)
{
    // Define the scalar type.
    // typedef float Scalar;
    // typedef double Scalar;
    // typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);

    // Define the functor.
    Scalar alpha = Scalar(1.);
    IntegrandExampleFunctor<Scalar> inFctr(alpha);

    // Define the integrator.
    Eigen::Integrator<Scalar> eigIntgtor(200);

    // Define a quadrature rule.
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod61;

    // Define the desired absolute and relative errors.
    Scalar desAbsErr = Scalar(0.);
    Scalar desRelErr = Eigen::NumTraits<Scalar>::epsilon() * Scalar(50.);

    // Integrate.
    Scalar result = eigIntgtor.quadratureAdaptive(inFctr, Scalar(0.), Scalar(1.), desAbsErr, desRelErr, quadratureRule);

    // Expected result.
    Scalar expected = Scalar(-4.);

    // Print output.
    int outputPrecision  = 18;
    std::cout << std::fixed;
    std::cout << "result          = " << std::setprecision(outputPrecision) << result << std::endl;
    std::cout << "exact result    = " << std::setprecision(outputPrecision) << expected << std::endl;
    std::cout << "actual error    = " << std::setprecision(outputPrecision) << (expected-result) << std::endl;

    return 0;
}
