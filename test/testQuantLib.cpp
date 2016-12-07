#include <NumericalIntegration.h>
#include <ql/math/functional.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/termstructures/volatility/abcd.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

/**
 * \param epsilon Relative machine precision.
 */
template <typename Scalar>
Scalar desiredRelativeError()
{
  return NumTraits<Scalar>::epsilon() * 50.;
}


void QuantLibSineIntegration(void)
{
    // Track the time required to complete the calculations.
    struct timeval timeStruct;
    gettimeofday(&timeStruct, NULL);
    long unsigned int processStartTime = timeStruct.tv_sec*1000000 + timeStruct.tv_usec;

    typedef QuantLib::Real Scalar;

    QuantLib::Size maxEvaluations = 1000;
    Scalar tolerance = desiredRelativeError<Scalar>();//1.e-6;

    std::cout<<"\nEvaluating QuantLib"<<std::endl;
    std::cout<<"Tolerance = "<<desiredRelativeError<Scalar>()<<"\n"<<std::endl;

    // 1) constant 0
    boost::function<Scalar (Scalar x)> fConst0;
    fConst0 = QuantLib::constant<Scalar,Scalar>(Scalar(0.));
    QuantLib::GaussKronrodAdaptive IConst0(tolerance,maxEvaluations);

    Scalar a(0.);
    Scalar b(1.);

    Scalar expected(0.);
    Scalar calculated = IConst0(fConst0,a,b);

    std::cout<<"constant 0,          |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 2) constant 1
    boost::function<Scalar (Scalar x)> fConst1;
    fConst1 = QuantLib::constant<Scalar,Scalar>(Scalar(1.));
    QuantLib::GaussKronrodAdaptive IConst1(tolerance,maxEvaluations);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(1.);
    calculated = IConst1(fConst1,a,b);

    std::cout<<"constant 1,          |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 3) identity
    boost::function<Scalar (Scalar x)> fIdentity;
    fIdentity = QuantLib::identity<Scalar>();
    QuantLib::GaussKronrodAdaptive IIdentity(tolerance,maxEvaluations);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(0.5);
    calculated = IIdentity(fIdentity,a,b);

    std::cout<<"identity,            |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 4) square
    boost::function<Scalar (Scalar x)> fSquare;
    fSquare = QuantLib::square<Scalar>();
    QuantLib::GaussKronrodAdaptive ISquare(tolerance,maxEvaluations);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(1.)/Scalar(3.);
    calculated = ISquare(fSquare,a,b);

    std::cout<<"square,              |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 5) sine
    boost::function<Scalar (Scalar x)> fSine;
    fSine = std::ptr_fun<Scalar,Scalar>(std::sin);
    QuantLib::GaussKronrodAdaptive ISine(tolerance,maxEvaluations);

    a = Scalar(0);
    b = Scalar(M_PI);

    expected = Scalar(2.);
    calculated = ISine(fSine,a,b);

    std::cout<<"sine,                |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 6) NormalDistribution
    boost::function<Scalar (Scalar x)> fNormalDistribution;
    fNormalDistribution = QuantLib::NormalDistribution();
    //maxEvaluations = 100000;
    QuantLib::GaussKronrodAdaptive INormalDistribution(tolerance,maxEvaluations);

    a = Scalar(-4);
    b = Scalar(4);

    expected = Scalar(1.);
    calculated = INormalDistribution(fNormalDistribution,a,b);

    std::cout<<"Normal Distribution, |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 7) Abcd2
    boost::function<Scalar (Scalar x)> fAbcd2;
    fAbcd2 = QuantLib::AbcdSquared(Scalar(0.07), Scalar(0.07),
        Scalar(0.5), Scalar(0.1), Scalar(8.0), Scalar(10.0));
    //maxEvaluations = 100000;
    QuantLib::GaussKronrodAdaptive IAbcd2(tolerance,maxEvaluations);

    a = Scalar(5);
    b = Scalar(6);

    expected = QuantLib::AbcdFunction(Scalar(0.07), Scalar(0.07), Scalar(0.5),
        Scalar(0.1)).covariance(Scalar(5.0),Scalar(6.0), Scalar(8.0), Scalar(10.0));
    calculated = IAbcd2(fAbcd2,a,b);

    std::cout<<"Abcd2,               |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    gettimeofday(&timeStruct, NULL);
    long unsigned int processFinishTime = timeStruct.tv_sec*1000000 + timeStruct.tv_usec;
    double totalTimeElapsed = (processFinishTime - processStartTime) / 1000000.;
    
    std::cout << "\n\tTotal Elapsed Time: " << totalTimeElapsed << std::endl;

}


void EigenSineIntegration()
{
    // Track the time required to complete the calculations.
    struct timeval timeStruct;
    gettimeofday(&timeStruct, NULL);
    long unsigned int processStartTime = timeStruct.tv_sec*1000000 + timeStruct.tv_usec;

    typedef QuantLib::Real Scalar;
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule =
        Eigen::Integrator<Scalar>::GaussKronrod15;
    const int maxSubintervals = 1000;

    std::cout<<"\nEvaluating Eigen "<<std::endl;
    std::cout<<"Tolerance = "<<desiredRelativeError<Scalar>()<<"\n"<<std::endl;

    typedef Eigen::Integrator<Scalar> IntegratorType;
    // 1) constant 0
    boost::function<Scalar (Scalar x)> fConst0;
    fConst0 = QuantLib::constant<Scalar,Scalar>(Scalar(0.));
    IntegratorType IConst0(maxSubintervals);

    Scalar a(0.);
    Scalar b(1.);

    Scalar expected(0.);
    Scalar calculated = IConst0.quadratureAdaptive(fConst0,a,b,Scalar(0.),
        desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"constant 0,          |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 2) constant 1
    boost::function<Scalar (Scalar x)> fConst1;
    fConst1 = QuantLib::constant<Scalar,Scalar>(Scalar(1.));
    IntegratorType IConst1(maxSubintervals);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(1.);
    calculated = IConst1.quadratureAdaptive(fConst1,a,b,Scalar(0.),
        desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"constant 1,          |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 3) identity
    boost::function<Scalar (Scalar x)> fIdentity;
    fIdentity = QuantLib::identity<Scalar>();
    IntegratorType IIdentity(maxSubintervals);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(0.5);
    calculated = IIdentity.quadratureAdaptive(fIdentity,a,b,Scalar(0.),
        desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"identity,            |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 4) square
    boost::function<Scalar (Scalar x)> fSquare;
    fSquare = QuantLib::square<Scalar>();
    IntegratorType ISquare(maxSubintervals);

    a = Scalar(0.);
    b = Scalar(1.);

    expected = Scalar(1.)/Scalar(3.);
    calculated = ISquare.quadratureAdaptive(fSquare,a,b,Scalar(0.),
        desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"square,              |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 5) sine
    boost::function<Scalar (Scalar x)> fSine;
    fSine = std::ptr_fun<Scalar,Scalar>(std::sin);
    IntegratorType ISine(maxSubintervals);

    a = Scalar(0);
    b = Scalar(M_PI);

    expected = Scalar(2.);
    calculated = ISine.quadratureAdaptive(fSine,a,b, Scalar(0.),
        desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"sine,                |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 6) NormalDistribution
    boost::function<Scalar (Scalar x)> fNormalDistribution;
    fNormalDistribution = QuantLib::NormalDistribution();
    //maxSubintervals = 100000;
    IntegratorType INormalDistribution(maxSubintervals);

    a = Scalar(-10);
    b = Scalar(10);

    expected = Scalar(1.);
    calculated = INormalDistribution.quadratureAdaptive(fNormalDistribution,a,b,
        Scalar(0.),desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"Normal Distribution, |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    // 7) Abcd2
    boost::function<Scalar (Scalar x)> fAbcd2;
    fAbcd2 = QuantLib::AbcdSquared(Scalar(0.07), Scalar(0.07),
        Scalar(0.5), Scalar(0.1), Scalar(8.0), Scalar(10.0));
    //maxEvaluations = 100000;
    IntegratorType IAbcd2(maxSubintervals);

    a = Scalar(5);
    b = Scalar(6);

    expected = QuantLib::AbcdFunction(Scalar(0.07), Scalar(0.07), Scalar(0.5),
        Scalar(0.1)).covariance(Scalar(5.0),Scalar(6.0), Scalar(8.0), Scalar(10.0));
    calculated = IAbcd2.quadratureAdaptive(fAbcd2,a,b,
        Scalar(0.),desiredRelativeError<Scalar>(),quadratureRule);

    std::cout<<"Abcd2,               |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

    gettimeofday(&timeStruct, NULL);
    long unsigned int processFinishTime = timeStruct.tv_sec*1000000 + timeStruct.tv_usec;
    double totalTimeElapsed = (processFinishTime - processStartTime) / 1000000.;
    
    std::cout << "\n\tTotal Elapsed Time: " << totalTimeElapsed << std::endl;
}

int main(void)
{
    QuantLibSineIntegration();
    EigenSineIntegration();
    return 0;
}
