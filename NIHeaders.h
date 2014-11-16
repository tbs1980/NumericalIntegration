#ifndef NI_NIHEADERS_H
#define NI_NIHEADERS_H

//#include <cassert>
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
T Pow(T base, T exponent)
{
    return pow(base,exponent);
}

float Pow(float base, float exponent)
{
    return std::pow(base,exponent);
}

double Pow(double base, double exponent)
{
    return std::pow(base,exponent);
}

long double Pow(long double base, long double exponent)
{
    return std::pow(base,exponent);
}

// Need to switch between gamma() from mpfrc++ and tgamma() from c++
template<typename T>
T Gamma(T x)
{
    return gamma(x);
}

double Gamma(double x)
{
    return tgamma(x);
}

float Gamma(float x)
{
    return tgamma(x);
}

long double Gamma(long double x)
{
    return tgamma(x);
}

#include "nodes_weights/kronrodLaurieGautschi.h"
#include "nodes_weights/kronrodPiessens.h"

#include "quadrature/QuadratureKronrod.h"
#include "quadrature/Integrator.h"

#endif //NI_NIHEADERS_H
