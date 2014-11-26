/* An example illustrating the use of numerical integration module in Eigen.
*/

// TODO change this header to something like NumericalInegration
// We may need to put all the code some common directory
#include <NIHeaders.h>
#include <iostream>
#include <iomanip>

/* We consider the example from
http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html

\int_0^1 x^{-1/2} log(x) dx = -4

Our integrator expects the user to provide a functor as shown below.

*/

template<typename Scalar>
class IntegrandExampleFunctor
{
public:
    IntegrandExampleFunctor(const Scalar alpha)
    :m_alpha(alpha)
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
    //define the sclar
    //typedef float Scalar;
    typedef double Scalar;
    //typedef long double Scalar;
    //typedef mpfr::mpreal Scalar;
    //Scalar::set_default_prec(256);

    //define the functor
    Scalar alpha=1.;
    IntegrandExampleFunctor<Scalar> inFctr(alpha);

    //define the integrator
    Eigen::Integrator<Scalar> eigIntgtor(200);

    //define a quadrature rule
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod61;

    //define the desired absolute and relative errors
    Scalar desAbsErr = Scalar(0.);
    Scalar desRelErr = Eigen::NumTraits<Scalar>::epsilon() * 50.;

    //integrate
    Scalar result = eigIntgtor.quadratureAdaptive(inFctr, Scalar(0.),Scalar(1.), desAbsErr, desRelErr, quadratureRule);

    //expected result
    Scalar expected = Scalar(-4.);

    //print output
    size_t outputPrecision  = 18;
    std::cout<<std::fixed;
    std::cout<<"result          = "<<std::setprecision(outputPrecision)<<result<<std::endl;
    std::cout<<"exact result    = "<<std::setprecision(outputPrecision)<<expected<<std::endl;
    std::cout<<"actual error    = "<<std::setprecision(outputPrecision)<<(expected-result)<<std::endl;


    return 0;
}
