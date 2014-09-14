#pragma once

using namespace Eigen;
using namespace std;

namespace Kronrod
{
void cSumInPlace(ArrayXd &arr);

Array<double,Dynamic,2> kronrod(const unsigned int nNodes, Array<double,Dynamic,2> ab);

Array<double,Dynamic,2> multiPrecisionKronrod(const unsigned int nNodes);

Array<double,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes);

Array<double,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes,double alpha);

Array<double,Dynamic,2> jacobiRecurrenceCoeff(const unsigned int nNodes, double alpha, double beta);

Array<double,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes);

Array<double,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes,double alpha);

Array<double,Dynamic,2> jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, double alpha, double beta);

Array<double,Dynamic,2> kronrodRecurrenceCoeff(const unsigned int nNodes, Array<double,Dynamic,2> ab0);
}

void Kronrod::cSumInPlace(ArrayXd &arr)
{
    size_t arrSize = arr.size();
    if(arrSize > 1)
    {
        for(size_t i = 1; i < arrSize ; i++)
        {
            arr(i) += arr(i - 1);
        }
    }
}

/**
*   multiPrecisionKronrod Arbitrary precision Kronrod abscissae & weights.
*
*   xGK = multiPrecisionKronrod(N) computes Kronrod points for (-1,1) with any required precision
*   using Multiprecision Computing Toolbox for MATLAB:
*   http://advanpix.com
*
*   Based on work of Dirk Laurie and Walter Gautschi.
*   Created by Pavel Holoborodko, November 7, 2011.
*/
Array<double,Dynamic,2> Kronrod::multiPrecisionKronrod(const unsigned int nNodes)
{
    Array<double,Dynamic,2> alphaBeta = jacobiRecurrenceCoeffZeroToOne(2 * nNodes);
    Array<double,Dynamic,2> g = kronrod(nNodes, alphaBeta);
    Array<double,Dynamic,2> xGK = ArrayX2d::Zero(2 * nNodes + 1, 2);

    xGK.col(0) = 2. * g.col(0) - 1.;
    xGK.col(1) = 2. * g.col(1);
    return xGK;
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
*   Supplied by Dirk Laurie, June 22, 1998.
*
*   Edited by Pavel Holoborodko, November 7, 2011:
*   Added arbitrary precision support using
*   Multiprecision Computing Toolbox for MATLAB http://advanpix.com
*   Ported to C++/Eigen by Grey Point Corporation August 2014.
*
*/
Array<double,Dynamic,2> Kronrod::kronrod(const unsigned int nNodes, Array<double,Dynamic,2> alphaBeta)
{
    Array<double,Dynamic,2> ab0 = kronrodRecurrenceCoeff(nNodes, alphaBeta);
    MatrixXd J = MatrixXd::Zero(2*nNodes + 1, 2*nNodes + 1);

    for(size_t i=0; i < 2*nNodes; ++i)
    {
        J(i,i) = ab0(i,0);
        J(i,i+1) = sqrt(ab0(i+1,1));
        J(i+1,i) = J(i, i+1);
    }

    J(2 * nNodes, 2 * nNodes) = ab0(2 * nNodes, 0);

    EigenSolver<MatrixXd> es(J);
    VectorXcd d = es.eigenvalues();
    MatrixXcd V = es.eigenvectors();

    //Bubble sort
    bool sorted = false;
    int i = 0;

    while(!sorted)
    {
        double di = d(i).real();
        double di1 = d(i+1).real();

        if(di1 < di)
        {
            complex<double> tmpD = d(i);
            d(i) = d(i+1);
            d(i+1) = tmpD;

            VectorXcd tmpV = V.col(i);
            V.col(i) = V.col(i+1);
            V.col(i+1) = tmpV;

            i = max(i-1, 0);
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

    ArrayXcd tempVec = V.col(V.cols() - 1);

    VectorXcd e = ab0(0,1) * (tempVec * tempVec);

    Array<double,Dynamic,2> xwGK = ArrayX2d::Zero(2*nNodes + 1, 2);
    xwGK.col(0) = d.real();
    xwGK.col(1) = e.real();


    return xwGK;
}

/**
*   R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
*
*   alphaBeta = r_jacobi(n, alpha, beta) generates the first n recurrence
*   coefficients for monic Jacobi polynomials with parameters
*   alpha and beta. These are orthogonal on [-1,1] relative to the
*   weight function w(t) = (1 - t)^a(1 + t)^b. The n alpha-coefficients
*   are stored in the first column, the n beta-coefficients in
*   the second column, of the nx2 array ab. The call ab=
*   R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
*   ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
*
*   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
*   Gautschi, 4-4-2002.
*/

Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes)
{
    return jacobiRecurrenceCoeff(nNodes, 0, 0);
}

Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes, double alpha)
{
    return jacobiRecurrenceCoeff(nNodes, alpha, alpha);
}

Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeff(const unsigned int nNodes, double alpha, double beta)
{
    double nu = (beta - alpha) / (alpha + beta + 2.);
    double mu = pow(2., alpha + beta + 1.) * tgamma(alpha + 1.) * tgamma(beta + 1.)
                / tgamma(alpha + beta + 2.);

    if (nNodes == 1)
    {
        Array<double,Dynamic,2> alphaBeta = ArrayXXd::Zero(1,2);
        alphaBeta(0,0) = nu;
        alphaBeta(0,1) = mu;

        return alphaBeta;
    }

    ArrayXd nAlphaBeta(nNodes);
    nAlphaBeta(0) = 1.;

    ArrayXd n(nNodes);

    for(size_t i = 0; i < nNodes; ++i)
    {
        n(i) = i;
    }

    nAlphaBeta = 2. * n + alpha + beta;

    ArrayXd A(nNodes);
    A = pow(beta,2.) - pow(alpha,2.) * ArrayXd::Ones(nNodes)
        / (nAlphaBeta * (nAlphaBeta + 2.));

    A(0) = nu;

    ArrayXd B(nNodes);
    B = 4. * (alpha + n) * (beta + n) * n * (alpha + beta + n)
        / ((nAlphaBeta * nAlphaBeta) * (nAlphaBeta + 1) * (nAlphaBeta - 1) );

    B(0) = mu;
    B(1) = 4. * (alpha + 1) * (beta + 1) / ( pow(alpha + beta + 2, 2) * (alpha + beta + 3) );
    //ab(rows,cols)

    Array<double,Dynamic,2> alphaBeta = ArrayXXd::Zero(nNodes,2);
    alphaBeta.col(0) = A;
    alphaBeta.col(1) = B;

    return alphaBeta;
}

/**
*   jacobiRecurrenceCoeffZeroToOne Recurrence coefficients for monic Jacobi polynomials
*   on [0,1].
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
*/
Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes)
{
    return jacobiRecurrenceCoeffZeroToOne(nNodes,0,0);
}

Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, double alpha)
{
    return jacobiRecurrenceCoeffZeroToOne(nNodes, alpha, alpha);
}

Array<double,Dynamic,2> Kronrod::jacobiRecurrenceCoeffZeroToOne(const unsigned int nNodes, double alpha, double beta)
{
    Array<double,Dynamic,2> coeffs = jacobiRecurrenceCoeff(nNodes, alpha, beta);
    Array<double,Dynamic,2> alphaBeta = ArrayXXd::Zero(nNodes,2);

    for (size_t i = 0; i < nNodes; i++)
    {
        alphaBeta(i, 0) = (1. + coeffs(i, 0))/2.;
    }

    alphaBeta(0, 1) = coeffs(0, 1) / pow(2., alpha + beta + 1.);

    for (size_t i = 1; i < nNodes; i++)
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
*   Supplied by Dirk Laurie, 6-22.1998
*
*   Edited by Pavel Holoborodko, November 7, 2011:
*   Added arbitrary precision support using
*   Multiprecision Computing Toolbox for MATLAB
*   http://advanpix.com
*/
Array<double,Dynamic,2> Kronrod::kronrodRecurrenceCoeff(const unsigned int nNodes, Array<double,Dynamic,2> ab0)
{
    ArrayXd alpha = ArrayXd::Zero(2 * nNodes + 1);
    ArrayXd beta = ArrayXd::Zero(2 * nNodes + 1);

    int i = 0;

    for (i = 0; i <= floor(3. * nNodes / 2.); ++i)
    {
        alpha(i) = ab0(i, 0);
    }

    for (i = 0; i <= ceil(3. * nNodes / 2.); ++i)
    {
        beta(i) = ab0(i, 1);
    }

    ArrayXd sig = ArrayXd::Zero(floor(nNodes / 2.) + 2);
    ArrayXd sigT = ArrayXd::Zero(floor(nNodes / 2.) + 2);

    sigT(1) = beta(nNodes + 1);

    // Loop 1
    for(size_t m = 0 ; m <= nNodes - 2; m++)
    {
        //k = floor((m+1)/2:-1:0 from matlab
        int kSize = floor((m + 1) / 2.) + 1;

        ArrayXd kIndexArr(kSize);
        ArrayXd lIndexArr(kSize);

        for(int jInc = kSize - 1; jInc >= 0 ; --jInc)
        {
            kIndexArr(kSize - 1 - jInc) = (double) jInc;
        }

        lIndexArr = m - kIndexArr;

        // Set up an array with coefficient-wise operations to be summed
        ArrayXd sumArray(kSize);

        for(i = 0; i < kIndexArr.size(); ++i)
        {
            sumArray(i) = (alpha(kIndexArr(i) + nNodes + 1) - alpha(lIndexArr(i))) * sigT(kIndexArr(i) + 1)
                          + beta(kIndexArr(i) + nNodes + 1) * sig(kIndexArr(i))
                          - beta(lIndexArr(i)) * sig(kIndexArr(i) + 1);
        }

        cSumInPlace(sumArray);

        //Assign sum array values to sigma array locations
        for(i = 0; i < kIndexArr.size(); ++i)
        {
            sig(kIndexArr(i) + 1) = sumArray(i);
        }

        ArrayXd tmpSig = sig;
        sig = sigT;
        sigT = tmpSig;
    }

    // Index array for sig shift
    int jVec1Size = floor(nNodes / 2.) + 1;
    ArrayXd jVec1(jVec1Size);

    for(i = jVec1Size - 1; i >= 0 ; --i)
    {
        jVec1(jVec1Size - 1 - i) = i;
    }

    // Shift sig contents
    ArrayXd newSig = sig;

    for(i = 0; i < jVec1Size ; ++i)
    {
        newSig(jVec1(i)+1) = sig(jVec1(i));
    }

    sig = newSig;

    // Loop 2
    for (size_t m = nNodes - 1; m <= 2 * nNodes - 3; ++m)
    {
        int kSize = floor((m-1.)/2.) - (m + 1 - nNodes) + 1;

        ArrayXd kIndexArr(kSize);
        ArrayXd lIndexArr(kSize);

        int index = 0;

        for (int kVal = m + 1 - nNodes; kVal <= floor((m-1.)/2.); kVal++)
        {
            kIndexArr(index) = (double) kVal;
            index++;
        }

        lIndexArr = m - kIndexArr;
        ArrayXd jIndexArr = nNodes - 1. - lIndexArr;

        // Set up an array with coefficient-wise operations to be summed
        ArrayXd sumArray(kSize);

        for (i = 0; i < kSize ; ++i)
        {
            sumArray(i) = -(alpha(kIndexArr(i)+nNodes + 1) - alpha(lIndexArr(i))) * sigT(jIndexArr(i) + 1)
                          - beta(kIndexArr(i)+nNodes+1)*sig(jIndexArr(i)+1)
                          + beta(lIndexArr(i)) * sig(jIndexArr(i)+2);
        }

        cSumInPlace(sumArray);

        //Assign sum array values to sigma array locations
        for (i = 0; i < kSize ; ++i)
        {
            sig(jIndexArr(i) + 1) = sumArray(i);
        }

        double jEnd = jIndexArr(kSize - 1);
        double kEnd = floor((m+1.)/2.);

        if (m % 2 == 0)
        {
            alpha(kEnd+nNodes + 1) = alpha(kEnd) + (sig(jEnd + 1) - beta(kEnd + nNodes + 1)
                                   * beta(kEnd + nNodes + 1) * sig(jEnd + 2)) / sigT(jEnd + 2);
        }
        else
        {
            beta(kEnd + nNodes + 1) = sig(jEnd + 1) / sig(jEnd + 2);
        }

        ArrayXd tmp = sig;
        sig = sigT;
        sigT = tmp;
    }

    alpha(2 * nNodes) = alpha(nNodes - 1) - beta(2 * nNodes) * sig(1) / sigT(1);

    Array<double,Dynamic, 2> alphaBeta = ArrayX2d::Zero(2 * nNodes + 1, 2);
    alphaBeta.col(0) = alpha.array();
    alphaBeta.col(1) = beta.array();

    return alphaBeta;
}
