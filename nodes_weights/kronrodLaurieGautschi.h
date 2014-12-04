/**
* The functions contained in this file calculate the Gauss-Kronrod nodes and weights
* using the Laurie/Gautschi method.
*/

#ifndef NI_KRONRODLAURIEGAUTSCHI_H
#define NI_KRONRODLAURIEGAUTSCHI_H

namespace Kronrod {

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
    template<typename _RealType>
    class LaurieGautschi
    {
    public:
        typedef _RealType RealType;
        typedef typename Eigen::Matrix<RealType,Eigen::Dynamic,1>  VectorType;
        typedef typename Eigen::Matrix<RealType,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
        typedef typename VectorType::Index IndexType;
        typedef typename Eigen::SelfAdjointEigenSolver<MatrixType> SelfAdjointEigenSolverType;

        /**
         * \brief Recurrence coefficients for monic Jacobi polynomials.
         *
         * This method generates the first N recurrence
         * coefficients for monic Jacobi polynomials with parameters
         * a and b. These are orthogonal on [-1,1] relative to the
         * weight function w(t)=(1-t)^a(1+t)^b. The N alpha-coefficients
         * are stored in \a a_out, the n beta-coefficients in \a b_out.
         * http://en.wikipedia.org/wiki/Jacobi_polynomials
         *
         * Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and  Matt Beall, September 2014
         *
         * \param[in] N Number of recurrence coefficients
         * \param[in] a alpha parameter of the Jacobi-polynomials
         * \param[in] b beta parameter of the Jacobi-polynomials
         * \param[in/out] a_out N alpha-coefficients
         * \param[in/out] b_out N beta-coefficients
         */
        static void r_jacobi(const IndexType N,const RealType a,const RealType b,
            VectorType & a_out, VectorType & b_out)
        {
            //TODO : make use the eigen assert facilities
            assert(a > RealType(-1));
            assert(b > RealType(-1));
            assert(a_out.rows()==b_out.rows());
            assert(a_out.rows() > 0);
            assert(N<=a_out.rows());

            using std::pow;
            a_out(0) = (b-a)/(a+b+RealType(2.));
            b_out(0) = pow(RealType(2.),(a+b+RealType(1.)))*Gamma(a+RealType(1.))*Gamma(b+RealType(1.))/Gamma(a+b+RealType(2.));

            for(IndexType n=1;n<N;++n)
            {
                RealType nab = RealType(2.)*n+a+b;
                a_out(n) = (b*b - a*a)/(nab*(nab+RealType(2.)));
                b_out(n) =  RealType(4.)*(n+a)*(n+b)*n*(n+a+b)/(nab*nab*(nab+RealType(1.))*(nab-RealType(1.)));
            }
        }

        /**
         * \brief Recurrence coefficients for monic Jacobi polynomials on [0,1].
         *
         * This method generates the first N recurrence
         * coefficients for monic Jacobi polynomials on [0,1] with
         * parameters a and b. These are orthogonal on [0,1] relative
         * to the weight function w(t)=(1-t)^a t^b. The N alpha-
         * coefficients are stored in \a a_out, the N beta-
         * coefficients in \a b_out.
         * http://en.wikipedia.org/wiki/Jacobi_polynomials
         *
         * Created by Dirk Laurie, 6-22-1998; edited by Walter Gautschi, 4-4-2002.
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and  Matt Beall, September 2014
         *
         * \param[in] N Number of recurrence coefficients
         * \param[in] a alpha parameter of the Jacobi-polynomials
         * \param[in] b beta parameter of the Jacobi-polynomials
         * \param[in/out] a_out N alpha-coefficients
         * \param[in/out] b_out N beta-coefficients
         */
        static void r_jacobi_01(const IndexType N,const RealType a,const RealType b,
            VectorType & a_out, VectorType & b_out)
        {
            //TODO : make use the eigen assert facilities
            assert(a > RealType(-1));
            assert(b > RealType(-1));
            assert(a_out.rows()==b_out.rows());
            assert(a_out.rows() > 0);
            assert(N<=a_out.rows());

            r_jacobi(N,a, b, a_out, b_out);

            for(IndexType n=0;n<N;++n)
            {
                a_out(n)  = (RealType(1.)+a_out(n))/RealType(2.);
            }

            using std::pow;
            b_out(0) = b_out(0)/pow(RealType(2),a+b+RealType(1.));

            for(IndexType n=1;n<N;++n)
            {
                b_out(n) = b_out(n)/RealType(4.);
            }

        }

        /**
         * \brief Jacobi-Kronrod matrix.
         *
         * This method produces the alpha- and beta-elements in
         * the Jacobi-Kronrod matrix of order 2N+1 for the weight
         * function (or measure) w. The input data for the weight
         * function w are the recurrence coefficients of the associated
         * orthogonal polynomials, which are stored in \a a_in and \a b_in .
         * At least ceil(3*N/2)+1 coefficients should be provided.
         * The 2N+1 alpha- and beta-elements are returned in \a a and \a b
         * respectively.
         *
         * Created by Dirk Laurie, 6-22.1998
         * Edited by Pavel Holoborodko, November 7, 2011
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder and
         * Matt Beall, September 2014
         *
         * \param[in] N Number of nodes
         * \param[in/out] a_in The recurrence coefficients of the associated orthogonal polynomials
         * \param[in/out] b_in The recurrence coefficients of the associated orthogonal polynomials
         * \param[in/out] a alpha-elements in the Jacobi-Kronrod matrix of order 2N+1
         * \param[in/out] b beta-elements in the Jacobi-Kronrod matrix of order 2N+1
         *
         */
        static void r_kronrod(const IndexType N,VectorType const & a_in, VectorType const & b_in,
            VectorType & a, VectorType & b)
        {
            //TODO : make use the eigen assert facilities
            assert(a_in.rows()==b_in.rows());
            assert(a_in.rows()>0);
            assert(a_in.rows() >= ceil(3*N/2)+1 );
            assert(a.rows()==2*N+1);
            assert(b.rows()==2*N+1);

            a=VectorType::Zero(2*N+1);
            b=VectorType::Zero(2*N+1);

            for(IndexType k=0;k<=floor(3*N/2)+1;++k)
            {
                a(k) = a_in(k);
            }

            for(IndexType k=0;k<=ceil(3*N/2)+1;++k)
            {
                b(k) = b_in(k);
            }

            VectorType s=VectorType::Zero(floor(N/2)+2);
            VectorType t=VectorType::Zero(floor(N/2)+2);

            t(1)=b(N+1);

            for(IndexType m=0;m<N-2+1;++m)
            {
                RealType u=0;
                for(IndexType k=floor((m+1)/2);k>=0;--k)
                {
                    IndexType l=m-k;
                    u = u + ( a(k+N+1)-a(l) )*t(k+1) + b(k+N+1)*s(k) - b(l)*s(k+1);
                    s(k+1) = u;
                }
                VectorType swap=s;
                s=t;
                t=swap;
            }

            for(IndexType j=floor(N/2);j>=0;--j)
            {
                s(j+1)=s(j);
            }

            for(IndexType m=N-1;m<2*N-3+1;++m)
            {
                IndexType k=m+1-N;
                IndexType j=0;
                RealType u=0;
                for(k=m+1-N;k<floor((m-1)/2)+1;++k)
                {
                    IndexType l=m-k;
                    j=N-1-l;
                    u = u - ( a(k+N+1)-a(l) )*t(j+1) - b(k+N+1)*s(j+1) + b(l)*s(j+2);
                    s(j+1) = u;
                }

                k=floor((m+1)/2);

                if(m % 2 == 0)
                {
                    a(k+N+1)=a(k)+(s(j+1)-b(k+N+1)*s(j+2))/t(j+2);
                }
                else
                {
                    b(k+N+1)=s(j+1)/s(j+2);
                }

                VectorType swap=s;
                s=t;
                t=swap;
            }

            a(2*N)=a(N-1)-b(2*N)*s(1)/t(1);
        }

        /**
         * \brief Gauss-Kronrod quadrature formula.
         *
         * This method generates the (2N+1)-point Gauss-Kronrod
         * quadrature rule for the weight function w encoded by the
         * recurrence matrix (a,b) of order [ceil(3*n/2)+1]x2 containing
         * in its first and second column respectively the alpha- and
         * beta-coefficients in the three-term recurrence relation
         * for w. The 2N+1 nodes, in increasing order, are output
         * into \a x, the corresponding weights into \a w .
         *
         * Created by Dirk Laurie, June 22, 1998.
         * Edited by Pavel Holoborodko, November 7, 2011:
         * Ported to C++/Eigen by Sreekumar Thaithara Balan, Mark Sauder
         * and  Matt Beall, September 2014
         *
         * \param[in] N Number of nodes
         * \param[in/out] a 2N alpha coefficients (input)
         * \param[in/out] b 2N beta coefficients (input)
         * \param[in/out] x 2N+1 nodes
         * \param[in/out] w 2N+1 weights corresponding to \a x
         */
        static void kronrod(const IndexType N,VectorType const & a, VectorType const & b,
            VectorType & x, VectorType & w)
        {
            //TODO : make use the eigen assert facilities
            assert(N>0);
            assert(a.rows()==2*N);
            assert(a.rows()==b.rows());
            assert(x.rows()==2*N+1);
            assert(x.rows()==w.rows());

            VectorType a0=VectorType::Zero(2*N+1);
            VectorType b0=VectorType::Zero(2*N+1);

            r_kronrod(N,a,b,a0,b0);

            //TODO : CHECK NEEDED LIKE THE ONE ON LINE 21 IN KRONROD.M
            // Do we have an approximately equal function in Eigen?
            assert( std::abs(b0.sum() - (RealType) (2*N+1)) > 1e-5 );

            MatrixType J=MatrixType::Zero(2*N+1,2*N+1);

            for(IndexType k=0;k<2*N;++k)
            {
                J(k,k)=a0(k);
                J(k,k+1)=sqrt(b0(k+1));
                J(k+1,k)=J(k,k+1);
            }

            J(2*N,2*N)=a0(2*N);

            //TODO : Is this assumption of positive definiteness correct?
            SelfAdjointEigenSolverType es(J);

            //TODO : make use the eigen assert facilities
            assert(es.info()==Eigen::Success);

            x=es.eigenvalues();
            MatrixType V=es.eigenvectors();

            w=b0(0)*(V.row(0).array()*V.row(0).array()).matrix();

        }

        /**
         * \brief Gauss quadrature rule.
         *
         * Given a weight function w encoded by (a,b) of the
         * first N recurrence coefficients for the associated orthogonal
         * polynomials, the first column of (a,b) containing the N alpha-
         * coefficients and the second column the N beta-coefficients,
         * the method generates the nodes and weights (x,w) of
         * the N-point Gauss quadrature rule for the weight function w.
         * The nodes, in increasing order, are stored in \a x ,
         * the N corresponding weights in \a w .
         *
         * \param[in] N Number of nodes
         * \param[in/out] a 2N alpha coefficients (input)
         * \param[in/out] b 2N beta coefficients (input)
         * \param[in/out] x 2N+1 nodes
         * \param[in/out] w 2N+1 weights corresponding to \a x
         */
        static void gauss(const IndexType N,VectorType const & a, VectorType const & b,
            VectorType & x, VectorType & w)
        {
            //TODO : make use the eigen assert facilities
            assert(N>0);
            assert(a.rows()==2*N);
            assert(a.rows()==b.rows());
            assert(x.rows()==N);
            assert(x.rows()==w.rows());

            MatrixType J=MatrixType::Zero(N,N);

            J(0,0)=a(0);
            for(IndexType n=1;n<N;++n)
            {
                J(n,n)=a(n);
                J(n,n-1)=sqrt(b(n));
                J(n-1,n)=J(n,n-1);
            }

            //TODO : Is this assumption of positive definiteness correct?
            SelfAdjointEigenSolverType es(J);

            //TODO : make use the eigen assert facilities
            assert(es.info()==Eigen::Success);

            x=es.eigenvalues();
            MatrixType V=es.eigenvectors();

            w=b(0)*(V.row(0).array()*V.row(0).array()).matrix();

        }

        /**
         * \brief Arbitrary precision Kronrod abscissae & weights.
         *
         * This method computes Kronrod points for (-1,1) with any required precision.
         * The result is a vector 2N+1 nodes (N points on either side of zero and zero)
         * and the corresponding weights.
         *
         * \param[in] N Number of nodes
         * \param[in/out] x Returns a vector of 2N+1 nodes
         * \param[in/out] w Returns a vector of weights corresponding to \a x
         */
        static void mpkronrod(const IndexType N,VectorType & x, VectorType & w)
        {
            //TODO : make use the eigen assert facilities
            assert(x.rows() ==  2*N+1);
            assert(w.rows() ==  2*N+1);
            assert(N>0);


            VectorType a=VectorType::Zero(2*N);
            VectorType b=VectorType::Zero(2*N);

            r_jacobi_01( 2*N, RealType(0), RealType(0), a, b);

            kronrod(N,a,b,x,w);

            for(IndexType i=0;i<x.rows();++i)
            {
                x(i) = RealType(2.)*x(i) - RealType(1.);
                w(i) = RealType(2.)*w(i);
            }
        }

        /**
         * \brief Arbitrary precision Gauss abscissae & weights.
         *
         * This method computes Kronrod points for (-1,1) with any required precision.
         * The result is a vector of N nodes the corresponding weights.
         *
         * \param[in] N Number of nodes
         * \param[in/out] x Returns a vector of 2N+1 nodes
         * \param[in/out] w Returns a vector of weights corresponding to \a x
         */
        static void mpgauss(const IndexType N,VectorType & x, VectorType & w)
        {
            //TODO : make use the eigen assert facilities
            assert(x.rows() ==  N);
            assert(w.rows() ==  N);
            assert(N>0);


            VectorType a=VectorType::Zero(2*N);
            VectorType b=VectorType::Zero(2*N);

            r_jacobi_01( 2*N, RealType(0), RealType(0), a, b);

            gauss(N,a,b,x,w);

            for(IndexType i=0;i<x.rows();++i)
            {
                x(i) = RealType(2.)*x(i) - RealType(1.);
                w(i) = RealType(2.)*w(i);
            }
        }

        static void computeAbscissaeAndWeights(unsigned int nNodes,
            Eigen::Array<RealType, Eigen::Dynamic, 1> &abscGaussKronrod,
            Eigen::Array<RealType, Eigen::Dynamic, 1> &weightGaussKronrod,
            Eigen::Array<RealType, Eigen::Dynamic, 1> &abscGauss,
            Eigen::Array<RealType, Eigen::Dynamic, 1> &weightGauss)
        {
            VectorType xGK = VectorType::Zero(2*nNodes+1);
            VectorType wGK = VectorType::Zero(2*nNodes+1);
            VectorType xG = VectorType::Zero(nNodes);
            VectorType wG = VectorType::Zero(nNodes);

            LaurieGautschi::mpkronrod(nNodes,xGK,wGK);
            LaurieGautschi::mpgauss(nNodes,xG,wG);

            unsigned int arraySize = nNodes + 1;

            abscGaussKronrod = Eigen::Array<RealType, Eigen::Dynamic, 1>::Zero(arraySize);
            weightGaussKronrod = Eigen::Array<RealType, Eigen::Dynamic, 1>::Zero(arraySize);
            abscGauss = Eigen::Array<RealType, Eigen::Dynamic, 1>::Zero(arraySize/2);
            weightGauss = Eigen::Array<RealType, Eigen::Dynamic, 1>::Zero(arraySize/2);

            using std::abs;
            for(unsigned int i=0;i<arraySize;++i)
            {
                abscGaussKronrod(i) = abs(xGK(i));
                weightGaussKronrod(i) = wGK(i);
            }

            abscGaussKronrod(arraySize-1) = RealType(0);

            for(unsigned int i=0;i<arraySize/2;++i)
            {
                abscGauss(i) = abs(xG(i));
                weightGauss(i) = wG(i);
            }
        }
    };
}

#endif //NI_KRONRODLAURIEGAUTSCHI_H
