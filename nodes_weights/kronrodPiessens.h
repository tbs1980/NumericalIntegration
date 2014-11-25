/**
* Piessens' implementation of the Kronrod calculations
*/

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cstdlib>
#include <iostream>

using namespace Eigen;
using namespace std;

namespace Kronrod
{
    template <typename Scalar>
    void kronrod(
        unsigned int nNodes, Array<Scalar, Dynamic, 1>& abscGaussKronrod,
        Array<Scalar, Dynamic, 1>& weightGaussKronrod, Array<Scalar,
        Dynamic, 1>& weightGauss);

    template <typename Scalar>
    void abscWeightKronrod(
        unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
        Array<Scalar, Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
        Scalar& weightGaussKronrod);

    template <typename Scalar>
    void abscWeightGauss(
        unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
        Array<Scalar,Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
        Scalar& weightGaussKronrod, Scalar& weightGauss);

    template <typename Scalar>
    Scalar machineEpsilon();
}

/**
* \brief kronrod adds n+1 points to an n-point Gaussian rule.
*
*    This function is a C++ implementation of original work by R. Piessens, et.al,
*    published in the journal Mathematics of Computation, Volume 28, Number 125,
*    January, 1974.  Where possible strucutre and variable naming convention
*    has been aligned to favor work of the QUADPACK Gaus-Kronrod integration
*    routines created by individuals of the same group.
*
*    This function calculates the abscissas and weights of the (2n+1)-point
*    Gauss-Kronrod quadrature formula which is obtained from the n-point
*    Gauss quadrature formula with the optimal addition of (n+1)-points.
*
*    The optimally added points are the Kronrod abscissae.  The
*    abscissas and weights for both the Gauss and Gauss Kronrod rules
*    are calculated for integration over the interval (-1, +1).
*
*    Because the quadrature formula is symmetric with respect to the origin,
*    only the positive abscissas are calculated.  Weights corresponding to the
*    symetric abscissae are equal.  Weights of weightGauss are calculated as well.
*
*    Work by Dr. John Burkhardt made note that the code published in Mathematics of
*    Computation omitted the definition of the second Chebyshev coefficient (chebCoeff2),
*    and Dr. Burkhardt's contributions are reflected here with permission.
*
*    The arrays abcsGaussKronrod, weightGaussKronrod and weightGauss contain the
*    positive abscissae in decreasing order, and the weights of each abscissa in
*    the Gauss-Kronrod and Gauss rules, respectively.
*
* Input Parameters:
* \param n, the order of the Gauss rule.
*
* Return Parameters:
* \return abscGaussKronrod[n+1] The Gauss-Kronrod abscissae.
* \return weightGaussKronrod[n+1] The weights for the Gauss-Kronrod rule.
* \return weightGauss[n+1] The weights for the Gauss rule.
*/
template <typename Scalar>
void Kronrod::kronrod(
    unsigned int nNodes, Array<Scalar, Dynamic, 1> &abscGaussKronrod,
    Array<Scalar, Dynamic, 1> &weightGaussKronrod, Array<Scalar, Dynamic, 1> &weightGauss)
{
    typedef Array<Scalar,Dynamic,1> ArrayXdType;
    unsigned int arraySize = nNodes + 1;
    abscGaussKronrod = ArrayXdType::Zero(arraySize);
    weightGaussKronrod = ArrayXdType::Zero(arraySize);
    weightGauss = ArrayXdType::Zero(arraySize / 2);

    Scalar aN = 0.0;
    Scalar d = 2.0;

    for (size_t i = 0; i < nNodes; ++i)
    {
        aN += 1.0;
        d *= aN / (aN + 0.5);
    }

    unsigned int m = (nNodes + 1) / 2;
    bool even = (nNodes == 2 * m);

    // aK is an index variable to account for only calculating positive abscissae
    Scalar aK = aN;

    // Calculation of the Chebyshev coefficients of the orthogonal polynomial.
    ArrayXdType tau(m);
    tau(0) = (aN + 2.0) / (2 * aN + 3.0);

    ArrayXdType betaCoeffs(m + 1);
    betaCoeffs(m - 1) = tau(0) - 1.0;

    for (size_t k = 1; k < m; ++k)
    {
        // This step accounts for both positive and negative abscissae
        aK += 2.0;

        tau(k) = ((aK - 1.0) * aK - aN * (aN + 1.0)) * (aK + 2.0) * tau(k-1) /
                 (aK * ((aK + 3.0) * (aK + 2.0) - aN * (aN + 1.0)));

        betaCoeffs(m-k-1) = tau(k);

        for (size_t i = 1; i <= k; ++i)
        {
            betaCoeffs(m - k - 1) = betaCoeffs(m - k + i - 1) * tau(i - 1) + betaCoeffs(m - k - 1);
        }
    }

    betaCoeffs(m) = 1.;

    // Calculation of approximate values for the abscissae as inital values
    // for the Newton-Raphson iterative solution.  These values are derived
    // from Pythagorean identities to the original code to more closely follow
    // the mathematics of the 1974 CoM paper.
    Scalar s1 = sin((M_PI / 2) / (2. * aN + 1.0));
    Scalar c1 = cos((M_PI / 2) / (2. * aN + 1.0));

    Scalar s2 = sin((M_PI) / (2. * aN + 1.0));
    Scalar c2 = cos((M_PI) / (2. * aN + 1.0));

    // Coefficient for Gauss and Kronrod abscissae and weights
    Scalar chebCoeff1 = 1.0 - 1.0 / (8.0 * aN * aN) + 1.0 / (8.0 * aN * aN * aN);
    Scalar chebCoeff2 = 2.0 / (Scalar(2. * nNodes + 1));

    for (size_t i = 1; i <= nNodes; ++i)
    {
        chebCoeff2 =  4.0 * chebCoeff2 * i / (nNodes + i);
    }

    Scalar abscK = chebCoeff1 * c1;
    Scalar temp = 0.;

    // Calculation of the K-th (Kronrod) abscissa and the corresponding weight.
    for (size_t k = 0; k < nNodes; ++k)
    {
        abscWeightKronrod(nNodes, m, even, chebCoeff2, betaCoeffs, abscK, weightGaussKronrod(k));
        abscGaussKronrod(k) = abscK;
        ++k;

        temp = c1;
        c1 = temp * c2 - s1 * s2;
        s1 = temp * s2 + s1 * c2;
        abscK = chebCoeff1 * c1;

        // Calculation of the k+1 (Gauss) abscissa and the corresponding weights.
        abscWeightGauss(nNodes, m, even, chebCoeff2, betaCoeffs, abscK, weightGaussKronrod(k),
                        weightGauss(k/2));
        abscGaussKronrod(k) = abscK;

        temp = c1;
        c1 = temp * c2 - s1 * s2;
        s1 = temp * s2 + s1 * c2;
        abscK = chebCoeff1 * c1;
    }

    // Add a Kronrod abscissa at the origin if n is even.
    if (even)
    {
        abscWeightKronrod(nNodes, m, even, chebCoeff2, betaCoeffs, abscK, weightGaussKronrod(nNodes));
    }

    // Set the abscissa value at the origin to zero and exit the function.
    abscGaussKronrod(nNodes) = 0.0;
    return;
}

/**
* \brief abscWeightKronrod calculates a Kronrod abscissa and weight.
*
*  Input Parameters:
* \param betaCoeffs[m+1] The Chebyshev coefficients.
* \param chebCoeff A value needed to compute weights.
* \param even A boolean variable that is TRUE if n is even.
* \param n The order of the Gauss rule.
* \param m The value of ( n + 1 ) / 2.
*
*  Input/output:
* \param abscGaussKronrod An estimate for the abscissa on input and the computed abscissa on output.
* \param weightGaussKronrod The Gauss-Kronrod weight.
*/
template <typename Scalar>
void Kronrod::abscWeightKronrod(
    unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
    Array<Scalar, Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
    Scalar& weightGaussKronrod)
{
    Scalar ai;

    Scalar b0 = 0.;
    Scalar b1 = 0.;
    Scalar b2 = 0.;

    Scalar d0 = 0.;
    Scalar d1 = 0.;
    Scalar d2 = 0.;

    Scalar delta = 1.;
    Scalar dif = 0.;

    Scalar f = 0.;
    Scalar fd = 0.;

    Scalar yy = 0.;

    int i = 0;

    size_t iter = 0;
    size_t iterationLimit = 50;

    // Iterative process for the computation of a Kronrod abscissa.
    while (fabs(delta) > machineEpsilon<Scalar>())
    {
        ++iter;

        b1 = 0.0;
        b2 = betaCoeffs(m);

        yy = 4.0 * abscGaussKronrod * abscGaussKronrod - 2.0;
        d1 = 0.0;

        if (even)
        {
            ai = m + m + 1;
            d2 = ai * betaCoeffs(m);
            dif = 2.0;
        }
        else
        {
            ai = m + 1;
            d2 = 0.0;
            dif = 1.0;
        }

        for (size_t k = 0; k < m; ++k)
        {
            ai -= dif;
            i = m - k - 1;
            b0 = b1;
            b1 = b2;
            d0 = d1;
            d1 = d2;
            b2 = yy * b1 - b0 + betaCoeffs(i);

            if (!even)
            {
                i += 1;
            }

            d2 = yy * d1 - d0 + ai * betaCoeffs(i);
        }

        if (even)
        {
            f = abscGaussKronrod * (b2 - b1);
            fd = d2 + d1;
        }
        else
        {
            f = 0.5 * (b2 - b0);
            fd = 4.0 * abscGaussKronrod * d2;
        }

        // Newton correction.
        delta = f / fd;
        abscGaussKronrod -= delta;

        if (abscGaussKronrod == Scalar(0.))
        {
            break;
        }

        // Identify non-convergence of the iterative solver after 50 iterations
        if (iter > iterationLimit)
        {
            cout << "Newton-Raphson iterative abscissae solver failed.";
            return;
        }
    }

    // Computation of the weight.
    d0 = 1.;
    d1 = abscGaussKronrod;
    ai = 0.;

    for (size_t k = 0; k < nNodes - 1; ++k)
    {
        ai = ai + 1.;
        d2 = ((ai + ai + 1.) * (abscGaussKronrod * d1) - (ai * d0)) / (ai + 1.);
        d0 = d1;
        d1 = d2;
    }

    weightGaussKronrod = chebCoeff / (fd * d2);
    return;
}

/**
* \brief abscWeightGauss calculates a Gaussian abscissa and two weights.
*
*   Input Parameters
* \param betaCoeffs[m+1] The Chebyshev coefficients.
* \param chebCoeff A value needed to compute weights.
* \param even A boolean variable that is TRUE if n is even.
* \param nNodes The order of the Gauss rule.
* \param m The value of ( nNodes + 1 ) / 2.
*
*   Input/Output
* \param abscGaussKronrod An estimate for the abscissa on input and the computed abscissa on output.
* \param weightGaussKronrod The Gauss-Kronrod weight.
* \param weightGauss The Gauss weight.
*/
template <typename Scalar>
void Kronrod::abscWeightGauss(
    unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
    Array<Scalar,Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
    Scalar& weightGaussKronrod, Scalar& weightGauss)
{
    Scalar ai = 0.;
    Scalar delta = 1.;

    Scalar p0 = 0.;
    Scalar p1 = 0.;
    Scalar p2 = 0.;

    Scalar pd0 = 0.;
    Scalar pd1 = 0.;
    Scalar pd2 = 0.;

    Scalar yy = 0.;

    size_t iter = 0;
    size_t iterationLimit = 50;

    //  Iterative process for the computation of a Gaussian abscissa.
    while (fabs(delta) > machineEpsilon<Scalar>())
    {
        ++iter;
        p0 = 1.;
        p1 = abscGaussKronrod;
        pd0 = 0.;
        pd1 = 1.;

        // If nNodes <= 1, initialize p2 and pd2 to avoid problems calculating delta.
        if (nNodes <= 1)
        {
            if (Kronrod::machineEpsilon<Scalar>() < fabs(abscGaussKronrod))
            {
                p2 = (3.0 * (abscGaussKronrod) * (abscGaussKronrod) - 1.0) / 2.0;
                pd2 = 3.0 * (abscGaussKronrod);
            }
            else
            {
                p2 = 3.0 * (abscGaussKronrod);
                pd2 = 3.0;
            }
        }

        ai = 0.0;

        for (size_t k = 0; k < nNodes - 1; ++k)
        {
            ai = ai + 1.0;
            p2 = ((ai + ai + 1.0) * abscGaussKronrod * p1 - ai * p0) / (ai + 1.0);
            pd2 = ((ai + ai + 1.0) * (p1 + abscGaussKronrod * pd1) - ai * pd0) / (ai + 1.0);
            p0 = p1;
            p1 = p2;
            pd0 = pd1;
            pd1 = pd2;
        }

        // Newton correction.
        delta = p2 / pd2;

        if (abscGaussKronrod == 0.0)
        {
            abscGaussKronrod -= delta;
            break;
        }
        else
        {
            abscGaussKronrod -= delta;
        }

        // Identify non-convergence of the iterative solver after 50 iterations
        if (iter > iterationLimit)
        {
            cout << "Newton-Raphson iterative abscissae solver failed.";
            return;
        }
    }

    // Computation of the Gauss weight.
    Scalar aN = nNodes;
    weightGauss = 2.0 / (aN * pd2 * p0);

    // Initialize
    p1 = 0.0;
    p2 = betaCoeffs(m);
    yy = 4.0 * (abscGaussKronrod) * (abscGaussKronrod) - 2.0;

    for (size_t k = 1; k <= m; ++k)
    {
        p0 = p1;
        p1 = p2;
        p2 = yy * p1 - p0 + betaCoeffs(m - k);
    }

    if (even)
    {
        weightGaussKronrod = weightGauss + chebCoeff / (abscGaussKronrod * pd2 * (p2 - p1));
    }
    else
    {
        weightGaussKronrod = weightGauss + (2.0 * chebCoeff) / (pd2 * (p2 - p0));
    }

    return;
}

/**
* \brief machineEpsilon returns the machine precision roundoff error.
* \returns Returns the machine precision round-off error for double precision
*/
template <typename Scalar>
Scalar Kronrod::machineEpsilon()
{
    // ISO standard machine precision values
    // half precision: 		2^-10	(9.76563e-04)
    // single precision: 	2^-23   (1.1920928955078125e-07)
    // double precision: 	2^-52   (2.22044604925031308e-16)
    // extended precision: 	2^-63	(1.08420217248550443e-19)
    // quad precision: 		2^-112	(1.9259299e-34)

    return Eigen::NumTraits<Scalar>::epsilon();
}
