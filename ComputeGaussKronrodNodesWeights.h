/**
* \file ComputeGaussKronrodNodesWeights.h
* The functions contained in this file calculate the Gauss-Kronrod nodes and weights
* using the Laurie/Gautschi method.
*/

#ifndef NI_KRONRODLAURIEGAUTSCHI_H
#define NI_KRONRODLAURIEGAUTSCHI_H

namespace Kronrod {

//----------------------------------Begin LaurieGautschi Class-----------------------------------//
    /**
     * \ingroup Quadrature_Module
     * \brief This class computes Kronrod abscissae & weights for arbitrary precision.
     *
     * Based on work of Dirk Laurie and Walter Gautschi.
     * D. P. Laurie (1997). Calculation of Gauss-Kronrod Quadrature Rules.
     * Mathematics of Computation, 66(219), 1133-1145
     * Created by Pavel Holoborodko, November 7, 2011.
     * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder, and Matt Beall September 2014
     *
     * \TODO Ensure only appropriates types are used for Scalar, e.g. prohibit integers.
     */
    template<typename _Scalar>
    class LaurieGautschi
    {
    public:
        typedef _Scalar Scalar;
        typedef typename Eigen::Matrix<Scalar,Eigen::Dynamic,1>  VectorType;
        typedef typename Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
        typedef typename VectorType::Index IndexType;
        typedef typename Eigen::SelfAdjointEigenSolver<MatrixType> SelfAdjointEigenSolverType;

        /**
         * \brief Recurrence coefficients for monic Jacobi polynomials.
         *
         * This method generates the first N recurrence
         * coefficients for monic Jacobi polynomials with parameters
         * alpha and beta. These are orthogonal on [-1,1] relative to the
         * weight function w(t)=(1-t)^a(1+t)^b. The N alpha-coefficients
         * are stored in \a alphaOut, the n beta-coefficients in \a betaOut.
         * http://en.wikipedia.org/wiki/Jacobi_polynomials
         *
         * Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and  Matt Beall, September 2014
         *
         * \param[in] N Number of recurrence coefficients
         * \param[in] alpha Alpha parameter of the Jacobi-polynomials
         * \param[in] beta Beta parameter of the Jacobi-polynomials
         * \param[in/out] alphaOut N alpha-coefficients
         * \param[in/out] betaOut N beta-coefficients
         */
        static void r_jacobi(const IndexType N,const Scalar alpha,const Scalar beta,
            VectorType& alphaOut, VectorType& betaOut)
        {
            //TODO : make use the eigen assert facilities
            assert(alpha > Scalar(-1));
            assert(beta > Scalar(-1));
            assert(alphaOut.rows() == betaOut.rows());
            assert(alphaOut.rows() > 0);
            assert(N <= alphaOut.rows());

            using std::pow;
            using std::tgamma;
            alphaOut(0) = (beta-alpha)/(alpha+beta+Scalar(2.));
            betaOut(0) = pow(Scalar(2.),(alpha+beta+Scalar(1.))) * tgamma(alpha+Scalar(1.)) * tgamma(beta+Scalar(1.)) / tgamma(alpha+beta+Scalar(2.));

            for(IndexType n=1;n<N;++n)
            {
                Scalar nAlphaBeta = Scalar(2.)*n+alpha+beta;
                alphaOut(n) = (beta*beta - alpha*alpha) / (nAlphaBeta*(nAlphaBeta+Scalar(2.)));
                betaOut(n) =  Scalar(4.)*(n+alpha)*(n+beta)*n*(n+alpha+beta) / (nAlphaBeta*nAlphaBeta*(nAlphaBeta+Scalar(1.))*(nAlphaBeta-Scalar(1.)));
            }
        }

        /**
         * \brief Recurrence coefficients for monic Jacobi polynomials on [0,1].
         *
         * This method generates the first N recurrence
         * coefficients for monic Jacobi polynomials on [0,1] with
         * parameters alpha and beta. These are orthogonal on [0,1] relative
         * to the weight function w(t)=(1-t)^alpha t^beta. The N alpha-
         * coefficients are stored in \a alphaOut, the N beta-
         * coefficients in \a betaOut.
         * http://en.wikipedia.org/wiki/Jacobi_polynomials
         *
         * Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and  Matt Beall, September 2014
         *
         * \param[in] N Number of recurrence coefficients
         * \param[in] alpha Alpha parameter of the Jacobi-polynomials
         * \param[in] beta Beta parameter of the Jacobi-polynomials
         * \param[in/out] alphaOut N alpha-coefficients
         * \param[in/out] betaOut N beta-coefficients
         */
        static void r_jacobi_01(const IndexType N,const Scalar alpha,const Scalar beta,
            VectorType& alphaOut, VectorType& betaOut)
        {
            //TODO : make use the eigen assert facilities
            assert(alpha > Scalar(-1));
            assert(beta > Scalar(-1));
            assert(alphaOut.rows() == betaOut.rows());
            assert(alphaOut.rows() > 0);
            assert(N <= alphaOut.rows());

            r_jacobi(N, alpha, beta, alphaOut, betaOut);

            for(IndexType n=0; n<N; ++n)
            {
                alphaOut(n) = (Scalar(1.)+alphaOut(n))/Scalar(2.);
            }

            using std::pow;
            betaOut(0) = betaOut(0)/pow(Scalar(2),alpha+beta+Scalar(1.));

            for(IndexType n=1; n<N; ++n)
            {
                betaOut(n) = betaOut(n)/Scalar(4.);
            }

        }

        /**
         * \brief Jacobi-Kronrod matrix.
         *
         * This method produces the alpha- and beta-elements in
         * the Jacobi-Kronrod matrix of order 2N+1 for the weight
         * function (or measure) w. The input data for the weight
         * function w are the recurrence coefficients of the associated
         * orthogonal polynomials, which are stored in \a alphaIn and \a betaIn .
         * At least ceil(3*N/2)+1 coefficients should be provided.
         * The 2N+1 alpha- and beta-elements are returned in \a alpha and \a beta
         * respectively.
         *
         * Created by Dirk Laurie, 6-22.1998
         * Edited by Pavel Holoborodko, November 7, 2011
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and  Matt Beall, September 2014
         *
         * \param[in] N Number of nodes
         * \param[in/out] alphaIn The recurrence coefficients of the associated orthogonal polynomials
         * \param[in/out] betaIn The recurrence coefficients of the associated orthogonal polynomials
         * \param[in/out] alpha Alpha-elements in the Jacobi-Kronrod matrix of order 2N+1
         * \param[in/out] beta Beta-elements in the Jacobi-Kronrod matrix of order 2N+1
         *
         */
        static void r_kronrod(const IndexType N,VectorType const& alphaIn, VectorType const& betaIn,
            VectorType& alpha, VectorType& beta)
        {
            //TODO : make use the eigen assert facilities
            assert(alphaIn.rows() == betaIn.rows());
            assert(alphaIn.rows() > 0);
            assert(alphaIn.rows() >= ceil(3*N/2)+1 );
            assert(alpha.rows() == 2*N+1);
            assert(beta.rows() == 2*N+1);

            alpha = VectorType::Zero(2*N+1);
            beta = VectorType::Zero(2*N+1);

            for(IndexType k=0; k<=floor(3*N/2)+1; ++k)
            {
                alpha(k) = alphaIn(k);
            }

            for(IndexType k=0; k<=ceil(3*N/2)+1; ++k)
            {
                beta(k) = betaIn(k);
            }

            VectorType sigma = VectorType::Zero(floor(N/2)+2);
            VectorType tempVector = VectorType::Zero(floor(N/2)+2);

            tempVector(1) = beta(N+1);

            for(IndexType m=0; m<N-2+1; ++m)
            {
                Scalar u = 0;
                for(IndexType k=floor((m+1)/2); k>=0;--k)
                {
                    IndexType l = m-k;
                    u = u + ( alpha(k+N+1)-alpha(l) )*tempVector(k+1) + beta(k+N+1)*sigma(k) - beta(l)*sigma(k+1);
                    sigma(k+1) = u;
                }

                VectorType swap = sigma;
                sigma = tempVector;
                tempVector = swap;
            }

            for(IndexType j=floor(N/2); j>=0; --j)
            {
                sigma(j+1) = sigma(j);
            }

            for(IndexType m = N-1; m<2*N-3+1; ++m)
            {
                IndexType k = m+1-N;
                IndexType j = 0;
                Scalar u = 0;
                for(k=m+1-N; k<floor((m-1)/2)+1; ++k)
                {
                    IndexType l = m-k;
                    j = N-1-l;
                    u = u - (alpha(k+N+1)-alpha(l))*tempVector(j+1) - beta(k+N+1)*sigma(j+1) + beta(l)*sigma(j+2);
                    sigma(j+1) = u;
                }

                k=floor((m+1)/2);

                if(m % 2 == 0)
                {
                    alpha(k+N+1) = alpha(k) + (sigma(j+1)-beta(k+N+1)*sigma(j+2)) / tempVector(j+2);
                }
                else
                {
                    beta(k+N+1) = sigma(j+1) / sigma(j+2);
                }

                VectorType swap = sigma;
                sigma = tempVector;
                tempVector = swap;
            }

            alpha(2*N) = alpha(N-1)-beta(2*N)*sigma(1)/tempVector(1);
        }

        /**
         * \brief Gauss-Kronrod quadrature formula.
         *
         * This method generates the (2N+1)-point Gauss-Kronrod
         * quadrature rule for the weight function w encoded by the
         * recurrence matrix (alpha,beta) of order [ceil(3*n/2)+1]x2 containing
         * in its first and second column respectively the alpha- and
         * beta-coefficients in the three-term recurrence relation
         * for w. The 2N+1 nodes, in increasing order, are output
         * into \a nodes, the corresponding weights into \a weights.
         *
         * Created by Dirk Laurie, June 22, 1998.
         * Edited by Pavel Holoborodko, November 7, 2011:
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder
         * and  Matt Beall, September 2014
         *
         * \param[in] N Number of nodes
         * \param[in/out] alpha 2N alpha coefficients (input)
         * \param[in/out] beta 2N beta coefficients (input)
         * \param[in/out] nodes 2N+1 nodes
         * \param[in/out] weights 2N+1 weights corresponding to \a nodes
         */
        static void kronrod(const IndexType N,VectorType const& alpha, VectorType const& beta,
            VectorType& nodes, VectorType& weights)
        {
            //TODO : make use the eigen assert facilities
            assert(N>0);
            assert(alpha.rows() == 2*N);
            assert(alpha.rows() == beta.rows());
            assert(nodes.rows() == 2*N+1);
            assert(nodes.rows() == weights.rows());

            VectorType alpha0 = VectorType::Zero(2*N+1);
            VectorType beta0 = VectorType::Zero(2*N+1);

            r_kronrod(N, alpha, beta, alpha0, beta0);

            // \TODO : CHECK NEEDED LIKE THE ONE ON LINE 21 IN KRONROD.M
            // Do we have an approximately equal function in Eigen?
            assert(std::abs(beta0.sum() - (Scalar) (2*N+1)) > 1e-5);

            MatrixType J = MatrixType::Zero(2*N+1,2*N+1);

            for(IndexType k=0; k<2*N; ++k)
            {
                J(k,k) = alpha0(k);
                J(k,k+1) = sqrt(beta0(k+1));
                J(k+1,k) = J(k,k+1);
            }

            J(2*N,2*N) = alpha0(2*N);

            //TODO : Is this assumption of positive definiteness correct?
            SelfAdjointEigenSolverType es(J);

            //TODO : make use the eigen assert facilities
            assert(es.info() == Eigen::Success);

            nodes = es.eigenvalues();
            MatrixType V = es.eigenvectors();

            weights = beta0(0)*(V.row(0).array()*V.row(0).array()).matrix();

        }

        /**
         * \brief Gauss quadrature rule.
         *
         * Given a weight function w encoded by (alpha,beta) of the
         * first N recurrence coefficients for the associated orthogonal
         * polynomials, the first column of (alpha,beta) containing the N alpha-
         * coefficients and the second column the N beta-coefficients,
         * the method generates the nodes and weights of
         * the N-point Gauss quadrature rule for the weight function.
         * The Nodes, in increasing order, are stored in \a nodes ,
         * the N corresponding weights are stored in \a weights .
         *
         * \param[in] N Number of nodes
         * \param[in/out] alpha 2N alpha coefficients (input)
         * \param[in/out] beta 2N beta coefficients (input)
         * \param[in/out] nodes 2N+1 nodes
         * \param[in/out] weights 2N+1 weights corresponding to \a nodes
         */
        static void gauss(const IndexType N,VectorType const& alpha, VectorType const& beta,
            VectorType& nodes, VectorType& weights)
        {
            //TODO : make use the eigen assert facilities
            assert(N > 0);
            assert(alpha.rows() == 2*N);
            assert(alpha.rows() == beta.rows());
            assert(nodes.rows() == N);
            assert(nodes.rows() == weights.rows());

            MatrixType J = MatrixType::Zero(N,N);

            J(0,0) = alpha(0);
            for(IndexType n=1; n<N; ++n)
            {
                J(n,n) = alpha(n);
                J(n,n-1) = sqrt(beta(n));
                J(n-1,n) = J(n,n-1);
            }

            //TODO : Is this assumption of positive definiteness correct?
            SelfAdjointEigenSolverType es(J);

            //TODO : make use the eigen assert facilities
            assert(es.info() == Eigen::Success);

            nodes = es.eigenvalues();
            MatrixType V = es.eigenvectors();

            weights = beta(0)*(V.row(0).array()*V.row(0).array()).matrix();

        }

        /**
         * \brief Arbitrary precision Kronrod abscissae & weights.
         *
         * This method computes Kronrod points for (-1,1) with any required precision.
         * The result is a vector 2N+1 nodes (N points on either side of zero and zero)
         * and the corresponding weights.
         *
         * \param[in] N Number of nodes
         * \param[in/out] nodes Returns a vector of 2N+1 nodes
         * \param[in/out] w Returns a vector of weights corresponding to \a nodes
         */
        static void mpkronrod(const IndexType N,VectorType& nodes, VectorType& weights)
        {
            //TODO : make use the eigen assert facilities
            assert(nodes.rows() ==  2*N+1);
            assert(weights.rows() ==  2*N+1);
            assert(N>0);

            VectorType alpha = VectorType::Zero(2*N);
            VectorType beta = VectorType::Zero(2*N);

            r_jacobi_01( 2*N, Scalar(0), Scalar(0), alpha, beta);

            kronrod(N,alpha,beta,nodes,weights);

            for(IndexType i=0; i<nodes.rows(); ++i)
            {
                nodes(i) = Scalar(2.)*nodes(i) - Scalar(1.);
                weights(i) = Scalar(2.)*weights(i);
            }
        }

        /**
         * \brief Arbitrary precision Gauss abscissae & weights.
         *
         * This method computes Kronrod points for (-1,1) with any required precision.
         * The result is a vector of N nodes the corresponding weights.
         *
         * \param[in] N Number of nodes
         * \param[in/out] nodes Returns a vector of 2N+1 nodes
         * \param[in/out] weights Returns a vector of weights corresponding to \a nodes
         */
        static void mpgauss(const IndexType N,VectorType& nodes, VectorType& weights)
        {
            //TODO : make use the eigen assert facilities
            assert(nodes.rows() == N);
            assert(weights.rows() == N);
            assert(N > 0);

            VectorType alpha = VectorType::Zero(2*N);
            VectorType beta = VectorType::Zero(2*N);

            r_jacobi_01(2*N, Scalar(0), Scalar(0), alpha, beta);

            gauss(N, alpha, beta, nodes, weights);

            for(IndexType i = 0;i<nodes.rows();++i)
            {
                nodes(i) = Scalar(2.)*nodes(i) - Scalar(1.);
                weights(i) = Scalar(2.)*weights(i);
            }
        }

        static void computeAbscissaeAndWeights(unsigned int nNodes,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& abscGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& weightGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& abscGauss,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& weightGauss)
        {
            VectorType xGK = VectorType::Zero(2*nNodes+1);
            VectorType wGK = VectorType::Zero(2*nNodes+1);
            VectorType xG = VectorType::Zero(nNodes);
            VectorType wG = VectorType::Zero(nNodes);

            LaurieGautschi::mpkronrod(nNodes,xGK,wGK);
            LaurieGautschi::mpgauss(nNodes,xG,wG);

            unsigned int arraySize = nNodes + 1;

            abscGaussKronrod = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(arraySize);
            weightGaussKronrod = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(arraySize);
            abscGauss = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(arraySize/2);
            weightGauss = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(arraySize/2);

            using std::abs;
            for(unsigned int i=0;i<arraySize;++i)
            {
                abscGaussKronrod(i) = abs(xGK(i));
                weightGaussKronrod(i) = wGK(i);
            }

            abscGaussKronrod(arraySize-1) = Scalar(0);

            for(unsigned int i=0;i<arraySize/2;++i)
            {
                abscGauss(i) = abs(xG(i));
                weightGauss(i) = wG(i);
            }
        }
    };
//-----------------------------------End LaurieGautschi Class------------------------------------//


//-------------------------------------Begin Piessens Class--------------------------------------//
    template<typename _Scalar>
    class Piessens
    {
    public:
        typedef _Scalar Scalar;

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
        *    Ported to C++/Eigen and templated for multiprecision by Mark Sauder,
        *    Sreekumar Thaithara Balan, Matt Beall, and R. Jeff Jenkins - September 2014.
        *
        * Input Parameters:
        * \param[in] n, the order of the Gauss rule.
        *
        * Return Parameters:
        * \param[in/out] abscGaussKronrod[n+1] The Gauss-Kronrod abscissae.
        * \param[in/out] weightGaussKronrod[n+1] The weights for the Gauss-Kronrod rule.
        * \param[in/out] weightGauss[n+1] The weights for the Gauss rule.
        */
        static void kronrod(
            unsigned int nNodes,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& abscGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& weightGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1>& weightGauss)
        {
            typedef Eigen::Array<Scalar,Eigen::Dynamic,1> ArrayXdType;
            unsigned int arraySize = nNodes + 1;
            abscGaussKronrod = ArrayXdType::Zero(arraySize);
            weightGaussKronrod = ArrayXdType::Zero(arraySize);
            weightGauss = ArrayXdType::Zero(arraySize / 2);

            Scalar aN(0.0);
            Scalar d(2.0);

            for (size_t i = 0; i < nNodes; ++i)
            {
                aN += Scalar(1.0);
                d *= aN / (aN + Scalar(0.5));
            }

            unsigned int m = (nNodes + 1) / 2;
            bool even = (nNodes == 2 * m);

            // aK is an index variable to account for calculating only the positive abscissae
            Scalar aK = aN;

            // Calculation of the Chebyshev coefficients of the orthogonal polynomial.
            ArrayXdType tau(m);
            tau(0) = (aN + Scalar(2.0)) / (Scalar(2) * aN + Scalar(3.0));

            ArrayXdType betaCoeffs(m + 1);
            betaCoeffs(m - 1) = tau(0) - Scalar(1.0);

            for (size_t k = 1; k < m; ++k)
            {
                // This step accounts for both positive and negative abscissae
                aK += Scalar(2.0);

                tau(k) = ((aK - Scalar(1.0)) * aK - aN * (aN + Scalar(1.0) )) * (aK + Scalar(2.0) ) * tau(k-1) /
                         (aK * ((aK + Scalar(3.0) ) * (aK + Scalar(2.0) ) - aN * (aN + Scalar(1.0) )));

                betaCoeffs(m-k-1) = tau(k);

                for (size_t i = 1; i <= k; ++i)
                {
                    betaCoeffs(m - k - 1) = betaCoeffs(m - k + i - 1) * tau(i - 1) + betaCoeffs(m - k - 1);
                }
            }

            betaCoeffs(m) = Scalar(1.);

            // Calculation of approximate values for the abscissae as inital values
            // for the Newton-Raphson iterative solution.  These values are derived
            // from Pythagorean identities to the original code to more closely follow
            // the mathematics of the 1974 CoM paper.

            // @TODO The usage of constant Pi with fixed precision needs to be changed to the following for multiprecision
            //RealScalar pi = NumTraits<RealScalar>::Pi();

            using std::sin;
            using std::cos;
            Scalar s1 = sin((M_PI / Scalar(2) ) / (Scalar(2.) * aN + Scalar(1.0) ));
            Scalar c1 = cos((M_PI / Scalar(2) ) / (Scalar(2.) * aN + Scalar(1.0) ));

            Scalar s2 = sin((M_PI) / (Scalar(2.) * aN + Scalar(1.0) ));
            Scalar c2 = cos((M_PI) / (Scalar(2.) * aN + Scalar(1.0) ));

            // Coefficient for Gauss and Kronrod abscissae and weights
            Scalar chebCoeff1 = Scalar(1.0) - Scalar(1.0) / (Scalar(8.0) * aN * aN) + Scalar(1.0) / (Scalar(8.0) * aN * aN * aN);
            Scalar chebCoeff2 = Scalar(2.0) / (Scalar(2. * nNodes + 1));

            for (size_t i = 1; i <= nNodes; ++i)
            {
                chebCoeff2 =  Scalar(4.0) * chebCoeff2 * i / (nNodes + i);
            }

            Scalar abscK = chebCoeff1 * c1;
            Scalar temp(0.);

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
            abscGaussKronrod(nNodes) = Scalar(0.0);
            return;
        }


        /**
        * \brief abscWeightKronrod calculates a Kronrod abscissa and weight.
        *
        *  Input Parameters:
        * \param[in] betaCoeffs[m+1] The Chebyshev coefficients.
        * \param[in] chebCoeff A value needed to compute weights.
        * \param[in] even A boolean variable that is TRUE if n is even.
        * \param[in] n The order of the Gauss rule.
        * \param[in] m The value of ( n + 1 ) / 2.
        *
        *  Input/output:
        * \param[in/out] abscGaussKronrod An estimate for the abscissa on input and the computed abscissa on output.
        * \param[in/out] weightGaussKronrod The Gauss-Kronrod weight.
        */
        static void abscWeightKronrod(
            unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
            Scalar& weightGaussKronrod)
        {
            Scalar ai;

            Scalar b0(0.);
            Scalar b1(0.);
            Scalar b2(0.);

            Scalar d0(0.);
            Scalar d1(0.);
            Scalar d2(0.);

            Scalar delta(1.);
            Scalar dif(0.);

            Scalar f(0.);
            Scalar fd(0.);

            Scalar yy(0.);

            int i = 0;

            size_t iter = 0;
            size_t iterationLimit = 50;

            // Iterative process for the computation of a Kronrod abscissa.
            using std::abs;
            while (abs(delta) > machineEpsilon())
            {
                ++iter;

                b1 = Scalar(0.0);
                b2 = betaCoeffs(m);

                yy = Scalar(4.0) * abscGaussKronrod * abscGaussKronrod - Scalar(2.0);
                d1 = Scalar(0.0);

                if (even)
                {
                    ai = Scalar(m + m + 1);
                    d2 = ai * betaCoeffs(m);
                    dif = Scalar(2.0);
                }
                else
                {
                    ai = Scalar(m + 1);
                    d2 = Scalar(0.0);
                    dif = Scalar(1.0);
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
                    f = Scalar(0.5) * (b2 - b0);
                    fd = Scalar(4.0) * abscGaussKronrod * d2;
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
                    std::cout << "@abscWeightKronrod Newton-Raphson iterative abscissae solver failed."<<std::endl;
                    return;
                }
            }

            // Computation of the weight.
            d0 = Scalar(1.);
            d1 = abscGaussKronrod;
            ai = Scalar(0.);

            for (size_t k = 0; k < nNodes - 1; ++k)
            {
                ai = ai + Scalar(1.);
                d2 = ((ai + ai + Scalar(1.)) * (abscGaussKronrod * d1) - (ai * d0)) / (ai + Scalar(1.));
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
        * \param[in] betaCoeffs[m+1] The Chebyshev coefficients.
        * \param[in] chebCoeff A value needed to compute weights.
        * \param[in] even A boolean variable that is TRUE if n is even.
        * \param[in] nNodes The order of the Gauss rule.
        * \param[in] m The value of ( nNodes + 1 ) / 2.
        *
        *   Input/Output
        * \param[in/out] abscGaussKronrod An estimate for the abscissa on input and the computed abscissa on output.
        * \param[in/out] weightGaussKronrod The Gauss-Kronrod weight.
        * \param[in/out] weightGauss The Gauss weight.
        */
        static void abscWeightGauss(
            unsigned int nNodes, unsigned int m, bool even, Scalar chebCoeff,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> betaCoeffs, Scalar& abscGaussKronrod,
            Scalar& weightGaussKronrod, Scalar& weightGauss)
        {
            Scalar ai(0.);
            Scalar delta(1.);

            Scalar p0(0.);
            Scalar p1(0.);
            Scalar p2(0.);

            Scalar pd0(0.);
            Scalar pd1(0.);
            Scalar pd2(0.);

            Scalar yy(0.);

            size_t iter = 0;
            size_t iterationLimit = 50;

            //  Iterative process for the computation of a Gaussian abscissa.
            using std::abs;
            while (abs(delta) > machineEpsilon())
            {
                ++iter;
                p0 = Scalar(1.);
                p1 = abscGaussKronrod;
                pd0 = Scalar(0.);
                pd1 = Scalar(1.);

                // If nNodes <= 1, initialize p2 and pd2 to avoid problems calculating delta.
                if (nNodes <= 1)
                {
                    using std::abs;
                    if (machineEpsilon() < abs(abscGaussKronrod))
                    {
                        p2 = (Scalar(3.0) * (abscGaussKronrod) * (abscGaussKronrod) - Scalar(1.0)) / Scalar(2.0);
                        pd2 = Scalar(3.0) * (abscGaussKronrod);
                    }
                    else
                    {
                        p2 = Scalar(3.0) * (abscGaussKronrod);
                        pd2 = Scalar(3.0);
                    }
                }

                ai = Scalar(0.0);

                for (size_t k = 0; k < nNodes - 1; ++k)
                {
                    ai = ai + Scalar(1.0);
                    p2 = ((ai + ai + Scalar(1.0) ) * abscGaussKronrod * p1 - ai * p0) / (ai + Scalar(1.0) );
                    pd2 = ((ai + ai + Scalar(1.0) ) * (p1 + abscGaussKronrod * pd1) - ai * pd0) / (ai + Scalar(1.0) );
                    p0 = p1;
                    p1 = p2;
                    pd0 = pd1;
                    pd1 = pd2;
                }

                // Newton correction.
                delta = p2 / pd2;

                if (abscGaussKronrod == Scalar(0.0))
                {
                    abscGaussKronrod -= delta;
                    break;
                }
                else
                {
                    abscGaussKronrod -= delta;
                }

                // Identify non-convergence of the iterative solver after iteration limit.
                if (iter > iterationLimit)
                {
                    std::cout << "@abscWeightGauss Newton-Raphson iterative abscissae solver failed."<<std::endl;;
                    return;
                }
            }

            // Computation of the Gauss weight.
            Scalar aN = nNodes;
            weightGauss = Scalar(2.0) / (aN * pd2 * p0);

            // Initialize
            p1 = Scalar(0.0);
            p2 = betaCoeffs(m);
            yy = Scalar(4.0) * (abscGaussKronrod) * (abscGaussKronrod) - Scalar(2.0);

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
                weightGaussKronrod = weightGauss + (Scalar(2.0) * chebCoeff) / (pd2 * (p2 - p0));
            }

            return;
        }

        /**
        * \brief machineEpsilon returns the machine precision roundoff error.
        * \returns Returns the machine precision round-off error for double precision
        */
        static Scalar machineEpsilon()
        {
            // ISO standard machine precision values
            // half precision:      2^-10   (9.76563e-04)
            // single precision:    2^-23   (1.1920928955078125e-07)
            // double precision:    2^-52   (2.22044604925031308e-16)
            // extended precision:  2^-63   (1.08420217248550443e-19)
            // quad precision:      2^-112  (1.9259299e-34)

            return Eigen::NumTraits<Scalar>::epsilon();
        }

        static void computeAbscissaeAndWeights(unsigned int nNodes,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> &abscGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> &weightGaussKronrod,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> abscGauss,
            Eigen::Array<Scalar, Eigen::Dynamic, 1> &weightGauss)
        {
            Piessens::kronrod(nNodes,abscGaussKronrod,weightGaussKronrod,weightGauss);
            abscGauss = Eigen::Array<Scalar, Eigen::Dynamic, 1>::Zero(abscGauss.rows());
        }
    };
//--------------------------------------End Piessens Class---------------------------------------//


//--------------------------------------Montegato Method---------------------------------------//

    /**
     * \brief A class for computign the Gauss/Kronrod weights.
     *
     * \tparam _Scalar floating point type
     */
    template<typename _Scalar>
    class Monegato
    {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ScalarArrayType;
        typedef typename ScalarArrayType::Index Index;

        static Scalar abs(Scalar x) {return (x < Scalar(0)) ? -(x) : x;}

        static Scalar legendre_err(const int n, const Scalar x, Scalar& err)
        {
            //using std::abs;
            if (n == 0)
            {
                err = Scalar(0);
                return Scalar(1);
            }
            else if (n == 1)
            {
                err = Scalar(0);
                return x;
            }

            // below Sree modified this to avoid -Wmaybe-uninitialized
            Scalar P0 = Scalar(1), P1 = x, P2=x;
            Scalar E0 = Eigen::NumTraits<Scalar>::epsilon();
            Scalar E1 = abs(x) * Eigen::NumTraits<Scalar>::epsilon();
            for (int k = 1; k < n; ++k)
            {
                P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
                err = ((2*k + 1) * abs(x) * E1 + k * E0) / (2*(k + 1));
                P0 = P1; P1 = P2;
                E0 = E1; E1 = err;
            }
            return P2;
        }

        static Scalar legendre_deriv(int const n, Scalar const x)
        {
            if (n == 0)
            {
                return Scalar(0);
            }
            else if (n == 1)
            {
                return Scalar(1);
            }

            Scalar P0 = Scalar(1), P1 = x, P2;
            // below Sree modified this to avoid -Wmaybe-uninitialized
            Scalar dP0 = Scalar(0), dP1 = Scalar(1), dP2=P1;
            for (int k = 1; k < n; ++k)
            {
                P2 = ((2*k + 1) * x * P1 - k * P0) / (k + Scalar(1));
                dP2 = (2*k + 1) * P1 + dP0;
                P0 = P1;
                P1 = P2;
                dP0 = dP1;
                dP1 = dP2;
            }
            return dP2;
        }

        static Scalar chebyshev_series_deriv(const Scalar x,const int n_,
            const ScalarArrayType& coefs)
        {
            Scalar d1(0), d2(0);
            Scalar y2 = 2 * x; // linear term for Clenshaw recursion

            for (int k = n_; k >= 2; --k)
            {
                Scalar temp = d1;
                d1 = y2 * d1 - d2 + k * coefs(k);
                d2 = temp;
            }

            return y2 * d1 - d2 + coefs(1);
        }

        static Scalar chebyshev_series(const Scalar x, const int n_,
            const ScalarArrayType& coefs, Scalar& err)
        {
            //using std::abs;
            Scalar d1(0), d2(0);
            Scalar absc = abs(coefs(0)); // final term for truncation error
            Scalar y2 = 2 * x; // linear term for Clenshaw recursion

            for (int k = n_; k >= 1; --k)
            {
                Scalar temp = d1;
                d1 = y2 * d1 - d2 + coefs(k);
                d2 = temp;
                absc += abs(coefs(k));
            }

            err = absc * Eigen::NumTraits<Scalar>::epsilon();
            return x * d1 - d2 + coefs(0)/2.;
        }


        static void legendre_zeros(const int m_,ScalarArrayType& zeros)
        {
            //using std::abs;
            ScalarArrayType temp = ScalarArrayType::Zero(m_+1);
            zeros(0) = Scalar(-1);
            zeros(1) = Scalar(1);
            Scalar delta, epsilon;

            for (int k = 1; k <= m_; ++k)
            {
                // loop to locate zeros of P_k interlacing z_0,...,z_k
                for (int j = 0; j < k; ++j)
                {
                    // Newton's method for P_k :
                    // initialize solver at midpoint of (z_j, z_{j+1})
                    delta = 1;
                    Scalar x_j = (zeros(j) + zeros(j+1)) / 2.;
                    Scalar P_k = legendre_err(k, x_j, epsilon);
                    while (abs(P_k) > epsilon &&
                        abs(delta) > Eigen::NumTraits<Scalar>::epsilon())
                    {
                        delta = P_k / legendre_deriv(k, x_j);
                        x_j -= delta;
                        P_k = legendre_err(k, x_j, epsilon);
                    }
                    temp(j) = x_j;
                }

                // copy roots tmp_0,...,tmp_{k-1} to z_1,...z_k:
                zeros(k+1) = zeros(k);
                for (int j = 0; j < k; ++j)
                {
                    zeros(j+1) = temp(j);
                }
            }

        }

        static void chebyshev_coefs(const int m_, ScalarArrayType& coefs)
        {
            int ell = (m_ + 1)/2;
            ScalarArrayType alpha = ScalarArrayType::Zero(ell+1);
            ScalarArrayType f = ScalarArrayType::Zero(ell+1);

            // Care must be exercised in initalizing the constants in the definitions.
            // Compilers interpret expressions like "(2*k + 1.0)/(k + 2.0)" as floating
            // point precision, before casting to Real.

            f(1) = Scalar(m_+1)/Scalar(2*m_ + 3);
            alpha(0) = Scalar(1); // coefficient of T_{m+1}
            alpha(1) = -f(1);

            for (int k = 2; k <= ell; ++k)
            {
                f(k) = f(k-1) * (2*k - 1) * (m_ + k) / (k * (2*m_ + 2*k + 1));
                alpha(k) = -f(k);
                for (int i = 1; i < k; ++i)
                {
                    alpha(k) -= f(i) * alpha(k-i);
                }
            }

            for (int k = 0; k <= ell; ++k)
            {
                coefs(m_ + 1 - 2*k) = alpha(k);
                if (m_  >= 2*k)
                {
                    coefs(m_ - 2*k) = Scalar(0);
                }
            }
        }

        static void gauss_kronrod_abscissae(const int n_,const int m_,
            const ScalarArrayType& zeros, const ScalarArrayType& coefs,
            ScalarArrayType& xgk_)
        {
            //
            // now from the function gauss_kronrod_abscissae
            //

            Scalar delta, epsilon;

            for (int k = 0; k < n_ / 2; ++k)
            {
                delta = 1;
                // Newton's method for E_{n+1} :
                Scalar x_k = (zeros(m_-k) + zeros(m_+1-k))/Scalar(2);
                Scalar E = chebyshev_series(x_k,n_,coefs, epsilon);
                while (abs(E) > epsilon &&
                    abs(delta) > Eigen::NumTraits<Scalar>::epsilon() )
                {
                    delta = E / chebyshev_series_deriv(x_k,n_,coefs);
                    x_k -= delta;
                    E = chebyshev_series(x_k,n_,coefs, epsilon);
                }
                xgk_(2*k) = x_k;
                // copy adjacent Legendre-zero into the array:
                if (2*k+1 < n_)
                {
                    xgk_(2*k+1) = zeros(m_-k);
                }
            }
        }

        static void gauss_kronrod_weights(const int& n_,const int m_,
            const ScalarArrayType& coefs, const ScalarArrayType& xgk_,
            ScalarArrayType& wg_, ScalarArrayType& wgk_ )
        {
            Scalar err;

            // Gauss-Legendre weights:
            for(int k=0; k < n_ /2; ++k)
            {
                Scalar x = xgk_(2*k + 1);
                wg_(k) = ( Scalar(-2)
                    / ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)) );
            }

            // The ratio of leading coefficients of P_n and T_{n+1} is computed
            // from the recursive formulae for the respective polynomials.
            Scalar F_m = Scalar(2) / Scalar(2*m_ + 1);
            for (int k = 1; k <= m_; ++k)
            {
                F_m *= (Scalar(2*k) / Scalar(2*k - 1));
            }

            // Gauss-Kronrod weights:

            for (int k = 0; k < n_; ++k)
            {
                Scalar x = xgk_(k);
                if (k % 2 == 0)
                {
                    wgk_(k) = F_m / (legendre_err(m_, x, err) * chebyshev_series_deriv(x,n_,coefs));
                }
                else
                {
                    wgk_(k) = (wg_(k/2) +
                        F_m / (legendre_deriv(m_, x) * chebyshev_series(x,n_,coefs, err)));
                }
            }

        }

        static void computeAbscissaeAndWeights(unsigned int m_,
            ScalarArrayType& xgk_,ScalarArrayType& wgk_,
            ScalarArrayType& xk_,ScalarArrayType& wg_)
        {
            //const Index m_ = 2*nNodes;
            const Index n_ = m_ + 1;

            std::cout<<"epsilon = "<<Eigen::NumTraits<Scalar>::epsilon()<<std::endl;

            xgk_ = ScalarArrayType::Zero(n_);//2*nNodes+1
            wgk_ = ScalarArrayType::Zero(n_);//2*nNodes+1
            xk_ = ScalarArrayType::Zero(n_/2);//2*nNodes
            wg_ = ScalarArrayType::Zero(n_/2);//2*nNodes

            // initialise the coefficients to zero
            ScalarArrayType coefs = ScalarArrayType::Zero(n_+1);
            ScalarArrayType zeros = ScalarArrayType::Zero(m_+2);

            legendre_zeros(m_, zeros);
            chebyshev_coefs(m_, coefs);
            gauss_kronrod_abscissae(n_, m_, zeros,coefs, xgk_);
            gauss_kronrod_weights(n_, m_,coefs, xgk_, wg_, wgk_ );
        }
    };
}

#endif //NI_KRONRODLAURIEGAUTSCHI_H
