#ifndef EIGEN_LAURIEGAUTSCHI_H
#define EIGEN_LAURIEGAUTSCHI_H

namespace Eigen {

    /**
     * \ingroup NumericalIntegration_Module
     *
     * \class LaurieGautschi
     *
     * \tparam _Scalar floating point type
     *
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
         * \param[in,out] alphaOut N alpha-coefficients
         * \param[in,out] betaOut N beta-coefficients
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
         * \param[in,out] alphaOut N alpha-coefficients
         * \param[in,out] betaOut N beta-coefficients
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
         * \param[in,out] alphaIn The recurrence coefficients of the associated orthogonal polynomials
         * \param[in,out] betaIn The recurrence coefficients of the associated orthogonal polynomials
         * \param[in,out] alpha Alpha-elements in the Jacobi-Kronrod matrix of order 2N+1
         * \param[in,out] beta Beta-elements in the Jacobi-Kronrod matrix of order 2N+1
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
         * \param[in,out] alpha 2N alpha coefficients (input)
         * \param[in,out] beta 2N beta coefficients (input)
         * \param[in,out] nodes 2N+1 nodes
         * \param[in,out] weights 2N+1 weights corresponding to \a nodes
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

            using std::sqrt;
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
         * \param[in,out] alpha 2N alpha coefficients (input)
         * \param[in,out] beta 2N beta coefficients (input)
         * \param[in,out] nodes 2N+1 nodes
         * \param[in,out] weights 2N+1 weights corresponding to \a nodes
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
            using std::sqrt;
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
         * \param[in,out] nodes Returns a vector of 2N+1 nodes
         * \param[in,out] w Returns a vector of weights corresponding to \a nodes
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
         * \param[in,out] nodes Returns a vector of 2N+1 nodes
         * \param[in,out] weights Returns a vector of weights corresponding to \a nodes
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

        /**
         * \brief A method ofr computing Gauss/Kronrod abscissae and weights
         * \param[in] nNodes Gauss-Legendre degree
         * \param[out] abscGaussKronrod Gauss/Kronrod abscissae
         * \param[out] weightGaussKronrod Gauss/Kronrod weights
         * \param[out] abscGauss Gauss abscissae
         * \param[out] weightGauss Gauss weights
         */
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

}//namespace Eigen

#endif //EIGEN_LAURIEGAUTSCHI_H
