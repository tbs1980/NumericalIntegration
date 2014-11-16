#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MPRealSupport>

template<typename T>
void test_pow(void)
{
    T base=0.1;
    T exponent=2.;

    T result = pow(base,exponent); //cannot use std::pow. mpreal compilation fails

    std::cout<<"result = "<<result<<std::endl;
}


namespace Eigen
{
    template<typename T>
    void test_pow(void)
    {
        T base=0.1;
        T exponent=2.;

        //T result = pow(base,exponent);//this will faill for float,double and long double
        T result = std::pow(base,exponent);//we need std::pow() if calling within Eigen

        std::cout<<"result = "<<result<<std::endl;
    }
}

int main(void)
{
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);

    test_pow<float>();
    test_pow<double>();
    test_pow<long double>();
    test_pow<Scalar>();


    Eigen::test_pow<float>();
    Eigen::test_pow<double>();
    Eigen::test_pow<long double>();
    //Eigen::test_pow<Scalar>(); /will not compile


    return 0;
}
