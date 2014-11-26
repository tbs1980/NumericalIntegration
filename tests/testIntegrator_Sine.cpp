#include <NIHeaders.h>

#include <iostream>
#include <fstream>
#include <iomanip>

template <typename Scalar>
Scalar desiredRelativeError()
{
  return Eigen::NumTraits<Scalar>::epsilon() * 50.;
}

template <typename Scalar>
typename Eigen::Integrator<Scalar>::QuadratureRule quadratureRules(const size_t i)
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

///////////////////////////// A simple example of sin(x) ///////////////////////
template<typename Scalar>
class IntegrandSineFunctor
{
public:
    Scalar operator()(const Scalar param) const
    {
        return sin(param);
    }
};


int test_sine(void)
{
    std::ofstream fout;
    fout.open("Sine_integration_test_output.txt");

    std::cout<<"Testing Int [0->Pi] sin(x) = 2"<<std::endl;

    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(117);

    typedef Eigen::Integrator<Scalar> IntegratorType;
    typedef IntegrandSineFunctor<Scalar> IntegrandSineFunctorType;

    //compute the nodes and weights on the fly
    QuadratureKronrod<Scalar>::computeNodesAndWeights();

    IntegratorType eigenIntegrator(256);
    IntegrandSineFunctorType integrandSineFunctor;

    bool success = true;
    const size_t numKeys = 12;

    for (size_t i = 0; i < numKeys; ++i)
    {
        success = true;

        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = quadratureRules<Scalar>(i);

        Scalar actual = eigenIntegrator.quadratureAdaptive(integrandSineFunctor, Scalar(0.), Scalar(M_PI),
            Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);

        Scalar expected = Scalar(2);

        if(fabs(expected - actual) > desiredRelativeError<Scalar>() * fabs(expected)
            or eigenIntegrator.errorCode() !=0)
        {
            fout << "\nrule " << i << "\n Abs(expected - actual) =" << fabs(expected - actual)
                      << "\n desiredRelativeError<Scalar>() * fabs(expected)= "
                      << desiredRelativeError<Scalar>() * fabs(expected) << std::endl;

            fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
            success = false;
        }
        else
        {
                fout << "\nrule " << i << "\n Abs(expected - actual) =" << fabs(expected - actual)
                          << "\n desiredRelativeError<Scalar>() * Abs(expected)= "
                          << desiredRelativeError<Scalar>() * fabs(expected) << std::endl;
                          
                fout << "Success!\n ";
        }

        if(success)    
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

    if (success)
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
    ret += test_sine();
    return EXIT_SUCCESS;
}
