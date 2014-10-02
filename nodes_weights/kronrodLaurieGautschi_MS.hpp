#ifndef EIGEN_QUADRATURE_KRONROD_MS_HPP
#define EIGEN_QUADRATURE_KRONROD_MS_HPP

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>

#include <iostream>

using namespace Eigen;

namespace Kronrod
{
    template <typename Scalar>
    Array<Scalar,Dynamic,2> multiPrecisionKronrod(const unsigned int nNodes);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> multiPrecisionGauss(const unsigned int nNodes);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes,Scalar alpha);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes, Scalar alpha, Scalar beta);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes,Scalar alpha);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, Scalar alpha, Scalar beta);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> gaussWeights(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> kronrod(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta);

    template <typename Scalar>
    Array<Scalar,Dynamic,2> kronrodRecurrenceCoeff(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta);

}

/**
*   multiPrecisionKronrod Arbitrary precision Kronrod abscissae & weights.
*
*   xwGK = multiPrecisionKronrod(N) computes Kronrod points for (-1,1) of any required precision
*
*   Based on work of Dirk Laurie and Walter Gautschi.
*   Created by Pavel Holoborodko, November 7, 2011.
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::multiPrecisionKronrod(const unsigned int nNodes)
{
    Array<Scalar,Dynamic,2> alphaBeta = jacobiRecurrenceCoeffZeroToOne<Scalar>(2 * nNodes);
    Array<Scalar,Dynamic,2> xwGK = kronrod(nNodes, alphaBeta);
    return xwGK;
}

/**
*   multiPrecisionKronrod Arbitrary precision Gauss abscissae & weights.
*
*   xwG = multiPrecisionKronrod(N) computes Kronrod points for (-1,1) of any required precision
*
*   Based on work of Dirk Laurie and Walter Gautschi.
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::multiPrecisionGauss(const unsigned int nNodes)
{
    Array<Scalar,Dynamic,2> alphaBeta = jacobiRecurrenceCoeffZeroToOne<Scalar>(2 * nNodes);
    Array<Scalar,Dynamic,2> xwG = gaussWeights(nNodes, alphaBeta);
    return xwG;
}


/**
*   gaussWeights Gauss quadrature formula.
*
*   Given a weight function w encoded by the Nx2 array alphaBeta of
*   the first nNodes recurrence coefficients for the associated orthogonal
*   polynomials, the first column of alphaBeta containing the N alpha-
*   coefficients and the second column the N beta-coefficients, the call
*   gaussWeights(nNodes,alphaBeta) generates the nodes and weights xwG
*   of the N-point Gauss quadrature rule for the weight function w.
*   The N nodes, in increasing order, are stored in the first column, the
*   N corresponding weights in the second column, of the Nx2 array xwG.
*   This method is often reffered to as the Golub-Welsch algorithm.
*
*   Based on work of Dirk Laurie and Walter Gautschi.
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::gaussWeights(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta)
{
    typedef Array<Scalar,Dynamic,1> ArrayXdType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorType;

    MatrixType jacobi = MatrixType::Zero(nNodes, nNodes);

    for (size_t n = 0; n < nNodes; ++n)
    {
        jacobi(n,n) = alphaBeta(n,0);
    }

    for (size_t n = 1; n < nNodes; ++n)
    {
      jacobi(n, n-1) = sqrt(alphaBeta(n,1));
      jacobi(n-1, n) = jacobi(n, n-1);
    }

    EigenSolver<MatrixType> eigenSol(jacobi);
    VectorType d = eigenSol.eigenvalues().real();
    MatrixType V = eigenSol.eigenvectors().real();

    // @TODO Is there a way to use std::sort()?
    //std::sort(d, d+d.size());
    //std::sort(V, V+V.size());

    //Insertion sort
    bool sorted = false;
    int i = 0;
    while(!sorted)
    {
        Scalar di = d(i);
        Scalar di1 = d(i+1);
        if(di1 < di)
        {
            Scalar tmpD = d(i);
            d(i) = d(i+1);
            d(i+1) = tmpD;
            VectorType tmpV = V.col(i);
            V.col(i) = V.col(i+1);
            V.col(i+1) = tmpV;
            i = std::max(i-1, 0);
            continue;
        }
        else
        {
            ++i;
            if(i == d.size() - 1)
            {
                sorted = true;
            }
        }
    }

    ArrayXdType tempV = V.row(0).array();
    ArrayXdType e = alphaBeta(0,1) * tempV * tempV;

    Array<Scalar,Dynamic,2> xwG = Array<Scalar,Dynamic,2>::Zero(nNodes,2);
    xwG.col(0) = d;
    xwG.col(1) = e;

    //Adjust the interval assumption of [0,1] to an interval of [-1, 1]
    xwG.col(0) = 2. * xwG.col(0) - 1.;
    xwG.col(1) = 2. * xwG.col(1);

    return xwG;
}

/**
*   kronrod Gauss-Kronrod quadrature formula.
*
*   xwGK = kronrod(n, alphaBeta) generates the (2n+1)-point Gauss-Kronrod
*   quadrature rule for the weight function w encoded by the
*   recurrence matrix alphaBeta of order [ceil(3*n/2)+1]x2 containing
*   in its first and second column respectively the alpha- and
*   beta-coefficients in the three-term recurrence relation
*   for w. The 2n+1 nodes, (abscissae), in increasing order, are output
*   into the first column. The corresponding weights are output to the
*   second column of the (2n+1)x2 array xwGK.
*
*   Created by Dirk Laurie, June 22, 1998.
*   Edited by Pavel Holoborodko, November 7, 2011:
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::kronrod(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta)
{
    typedef Array<Scalar,Dynamic,1> ArrayXdType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorType;

    Array<Scalar,Dynamic,2> alphaBeta0 = kronrodRecurrenceCoeff(nNodes, alphaBeta);
    MatrixType jacobi = MatrixType::Zero(2*nNodes + 1, 2*nNodes + 1);

    for(size_t k = 0; k < 2 * nNodes; ++k)
    {
        jacobi(k,k) = alphaBeta0(k,0);
        jacobi(k,k + 1) = sqrt(alphaBeta0(k + 1, 1));
        jacobi(k + 1, k) = jacobi(k, k + 1);
    }

    jacobi(2 * nNodes, 2 * nNodes) = alphaBeta0(2 * nNodes, 0);

    EigenSolver<MatrixType> eigenSol(jacobi);
    VectorType d = eigenSol.eigenvalues().real();
    MatrixType V = eigenSol.eigenvectors().real();

    // @TODO Is there a way to use std::sort()?
    //std::sort(d, d+d.size());
    //std::sort(V, V+V.size());

    //Insertion sort
    bool sorted = false;
    int i = 0;
    while(!sorted)
    {
        Scalar di = d(i);
        Scalar di1 = d(i+1);
        if(di1 < di)
        {
            Scalar tmpD = d(i);
            d(i) = d(i+1);
            d(i+1) = tmpD;
            VectorType tmpV = V.col(i);
            V.col(i) = V.col(i+1);
            V.col(i+1) = tmpV;
            i = std::max(i-1, 0);
            continue;
        }
        else
        {
            ++i;
            if(i == d.size() - 1)
            {
                sorted = true;
            }
        }
    }

    ArrayXdType tempV = V.row(0).array();
    ArrayXdType e = alphaBeta0(0,1) * tempV * tempV;

    Array<Scalar,Dynamic,2> xwGK = Array<Scalar,Dynamic,2>::Zero(2*nNodes + 1, 2);
    xwGK.col(0) = d.real();
    xwGK.col(1) = e;

    //Adjust the interval assumption of [0,1] to an interval of [-1, 1]
    xwGK.col(0) = 2. * xwGK.col(0) - 1.;
    xwGK.col(1) = 2. * xwGK.col(1);

    return xwGK;
}

/**
*   jacobiRecurrenceCoeff Recurrence coefficients for monic Jacobi polynomials.
*
*   alphaBeta = jacobiRecurrenceCoeff(n, alpha, beta) generates the first n
*   recurrence coefficients for monic Jacobi polynomials with parameters
*   alpha and beta. These are orthogonal on [-1,1] relative to the
*   weight function w(t) = (1 - t)^a(1 + t)^b. The n alpha-coefficients
*   are stored in the first column, the n beta-coefficients in the second column,
*   of the nx2 array ab. The call ab = jacobiRecurrenceCoeff(n,a)
*   is the same as ab = jacobiRecurrenceCoeff(n,a,a) and
*   ab = jacobiRecurrenceCoeff(n) is the same as ab = jacobiRecurrenceCoeff(n,0,0).
*
*   Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes)
{
    return jacobiRecurrenceCoeff(nNodes,(Scalar) 0,(Scalar) 0);
}

template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes, Scalar alpha)
{
    return jacobiRecurrenceCoeff(nNodes, alpha, alpha);
}

template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes, Scalar alpha, Scalar beta)
{
    typedef Array<Scalar,Dynamic,1> ArrayXdType;
    Scalar nu = (beta - alpha) / (alpha + beta + 2.);
    Scalar mu = pow(2., alpha + beta + 1.) * Gamma(alpha + 1.) * Gamma(beta + 1.)
                / Gamma(alpha + beta + 2.);

    if (nNodes == 1)
    {
        Array<Scalar,Dynamic,2> alphaBeta = Array<Scalar,Dynamic,2>::Zero(1,2);
        alphaBeta(0,0) = nu;
        alphaBeta(0,1) = mu;
        return alphaBeta;
    }

    Scalar  nAlphaBeta;
    ArrayXdType A(nNodes);
    ArrayXdType B(nNodes);

    for(size_t n = 1; n < nNodes; ++n)
    {
        nAlphaBeta = 2. * n + alpha + beta;
        A(n) = (pow(beta, 2.) - pow(alpha, 2.)) / (nAlphaBeta * (nAlphaBeta + 2.));

        B(n) = 4. * (alpha + n) * (beta + n) * n * (alpha + beta + n)
            / ((nAlphaBeta * nAlphaBeta) * (nAlphaBeta + 1.) * (nAlphaBeta - 1.));
    }

    B(0) = mu;
    B(1) = 4. * (alpha + 1.) * (beta + 1.) / (pow(alpha + beta + 2., 2.) * (alpha + beta + 3.));

    Array<Scalar,Dynamic,2> alphaBeta = Array<Scalar,Dynamic,2>::Zero(nNodes,2);
    alphaBeta.col(0) = A;
    alphaBeta.col(1) = B;

    return alphaBeta;
}

/**
*   jacobiRecurrenceCoeffZeroToOne Recurrence coefficients for monic Jacobi polynomials on [0,1].
*
*   alphaBeta = jacobiRecurrenceCoeffZeroToOne(n, a, b) generates the first n recurrence
*   coefficients for monic Jacobi polynomials on [0,1] with parameters alpha and beta.
*
*   These are orthogonal on [0,1] relative to the weight function w(t) = (1-t)^a t^b.
*   The n alpha-coefficients are stored in the first column, the n beta-coefficients
*   in the second column, of the nx2 array alphaBeta. The call
*   alphaBeta = jacobiRecurrenceCoeffZeroToOne(n, alpha) is the same as the call
*   alphaBeta = jacobiRecurrenceCoeffZeroToOne(n, alpha, alpha), and the call
*   alphaBeta = jacobiRecurrenceCoeffZeroToOne(n) is the same as the call
*   alphaBeta = jacobiRecurrenceCoeffZeroToOne(n, 0, 0).
*
*   Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes)
{
    return jacobiRecurrenceCoeffZeroToOne(nNodes,(Scalar)0,(Scalar)0);
}

template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, Scalar alpha)
{
    return jacobiRecurrenceCoeffZeroToOne(nNodes, alpha, alpha);
}

template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, Scalar alpha, Scalar beta)
{
    Array<Scalar,Dynamic,2> coeffs = jacobiRecurrenceCoeff(nNodes, alpha, beta);
    Array<Scalar,Dynamic,2> alphaBeta = Array<Scalar,Dynamic,2>::Zero(nNodes, 2);

    for (size_t i = 0; i < nNodes; ++i)
    {
        alphaBeta(i, 0) = (1. + coeffs(i, 0)) / 2.;
    }

    alphaBeta(0, 1) = coeffs(0, 1) / pow(2., alpha + beta + 1);

    for (size_t i = 1; i < nNodes; ++i)
    {
        alphaBeta(i, 1) = coeffs(i, 1) / 4.;
    }

    return alphaBeta;
}

/**
*   kronrodRecurrenceCoeff Jacobi-Kronrod matrix.
*
*   alphaBeta = kronrodRecurrenceCoeff(n,ab0) produces the alpha- and beta-elements in
*   the Jacobi-Kronrod matrix of order 2n+1 for the weight
*   function (or measure) w. The input data for the weight
*   function w are the recurrence coefficients of the associated
*   orthogonal polynomials, which are stored in the array ab0 of
*   dimension [ceil(3*n/2)+1]x2. The alpha-elements are stored in
*   the first column, the beta-elements in the second column, of
*   the (2*n+1)x2 array ab.
*
*   Created by Dirk Laurie, 6-22.1998
*   Edited by Pavel Holoborodko, November 7, 2011:
*   Ported to C++/Eigen Grey Point Corporation September 2014
*/
template <typename Scalar>
Array<Scalar,Dynamic,2> Kronrod::kronrodRecurrenceCoeff(const unsigned int nNodes, Array<Scalar,Dynamic,2> alphaBeta)
{
    typedef Array<Scalar,Dynamic,1> ArrayXdType;
    ArrayXdType alpha = ArrayXdType::Zero(2 * nNodes + 1);
    ArrayXdType beta = ArrayXdType::Zero(2 * nNodes + 1);

    Scalar temp = 0.;

    int j = 0;
    int k = 0;

    unsigned int m = 0;

    for (k = 0; k <= floor(3 * nNodes / 2 + 1); ++k)
    {
        alpha(k) = alphaBeta(k, 0);
    }

    for (k = 0; k <= ceil(3 * nNodes / 2 + 1); ++k)
    {
        beta(k) = alphaBeta(k, 1);
    }

    ArrayXdType sig = ArrayXdType::Zero(floor(nNodes / 2) + 2);
    ArrayXdType sigT = ArrayXdType::Zero(floor(nNodes / 2) + 2);

    sigT(1) = beta(nNodes + 1);

    // Loop 1
    for(m = 0; m < nNodes - 1; ++m)
    {
        temp = 0.;
        for(k = (m + 1) / 2; k >= 0; --k)
        {
            temp += (alpha(k + nNodes + 1) - alpha(m - k)) * sigT(k + 1)
                    + beta(k + nNodes + 1) * sig(k) - beta(m - k) * sig(k+1);
            sig(k+1) = temp;
        }

        ArrayXdType tempArray = sig;
        sig = sigT;
        sigT = tempArray;
    }

    for(j = floor(nNodes / 2); j >= 0; --j)
    {
        sig(j + 1) = sig(j);
    }

    // Loop 2
    for (m = nNodes - 1; m <= 2 * nNodes - 3; ++m)
    {
        temp = 0.;
        for(k = m - nNodes + 1; k <= floor((m - 1) / 2); ++k)
        {
            j = nNodes - m + k;
            temp = temp - (alpha(k + nNodes + 1) - alpha(m - k)) * sigT(j)
                   - beta(k + nNodes + 1) * sig(j)
                   + beta(m - k) * sig(j+1);
            sig(j) = temp;
        }

        if(m % 2 == 0)
        {
            alpha(k + nNodes + 1) = alpha(k) + (sig(j) - beta(k + nNodes + 1)
                                    * sig(j + 1)) / sigT(j + 1);
        }
        else
        {
            beta(k + nNodes + 1) = sig(j) / sig(j+1);
        }

        ArrayXdType tempArray = sig;
        sig = sigT;
        sigT = tempArray;
    }

    alpha(2 * nNodes) = alpha(nNodes - 1) - beta(2 * nNodes) * sig(1) / sigT(1);

    Array<Scalar,Dynamic, 2> alphaBeta0 = Array<Scalar,Dynamic,2>::Zero(2 * nNodes + 1, 2);
    alphaBeta0.col(0) = alpha.array();
    alphaBeta0.col(1) = beta.array();

    return alphaBeta0;
}

#endif //EIGEN_QUADRATURE_KRONROD_MS_HPP
