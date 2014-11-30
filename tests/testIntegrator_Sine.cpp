#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace Eigen;



#define NI_M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049

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
        return Sin(param);
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

        Scalar actual = eigenIntegrator.quadratureAdaptive(integrandSineFunctor, Scalar(0.), Scalar(NI_M_PI),
            Scalar(0.), desiredRelativeError<Scalar>(), quadratureRule);

        Scalar expected = Scalar(2);

        if(Abs(expected - actual) > desiredRelativeError<Scalar>() * Abs(expected)
            or eigenIntegrator.errorCode() !=0)
        {
            fout << "\nrule " << i << "\n Abs(expected - actual) =" << Abs(expected - actual)
                      << "\n desiredRelativeError<Scalar>() * Abs(expected)= "
                      << desiredRelativeError<Scalar>() * Abs(expected) << std::endl;

            fout << "errorCode = " << eigenIntegrator.errorCode() << std::endl;
            success = false;
        }
        else
        {
                fout << "\nrule " << i << "\n Abs(expected - actual) =" << Abs(expected - actual)
                          << "\n desiredRelativeError<Scalar>() * Abs(expected)= "
                          << desiredRelativeError<Scalar>() * Abs(expected) << std::endl;

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
