/* An example illustrating the use of numerical integration module in Eigen.
*/

// TODO change this header to something like NumericalInegration
// We may need to put all the codes some common directory
#include <NIHeaders.hpp>
#include <iostream>

/* We consider the example from
http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html#Numerical-integration-examples

\int_0^1 x^{-1/2} log(x) dx = -4

Our integrator expects the user to provide a functor as shown below.

*/

template<typename Scalar>
class IntegrandExampleFunctor
{
public:
    IntegrandExampleFunctor(const Scalar alpha)
    :m_alpha(alpha)

    Scalar operator()(const Scalar param) const
    {
        log(m_alpha*x) / sqrt(x);
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
    typedef double Scalar;

    //define the functor
    Scalar alpha=1.;
    IntegrandExampleFunctor<Scalar> itfctr(alpha);

    //define a quadrature rule
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod61;

    return 0;
}
