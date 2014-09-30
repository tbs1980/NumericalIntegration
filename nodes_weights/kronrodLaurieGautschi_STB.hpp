#ifndef NI_KRONRADLAURIEGAUTSCHI_STB_HPP
#define NI_KRONRADLAURIEGAUTSCHI_STB_HPP

namespace Kronrod {

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

    template<typename _RealType>
    class LaurieGautschi
    {
    public:
        typedef _RealType RealType;
        typedef typename Eigen::Matrix<RealType,Eigen::Dynamic,1>  VectorType;
        typedef typename Eigen::Matrix<RealType,Eigen::Dynamic, Eigen::Dynamic> MatrixType;
        typedef typename VectorType::Index IndexType;
        typedef typename Eigen::SelfAdjointEigenSolver<MatrixType> SelfAdjointEigenSolverType;

        static void r_jacobi(const IndexType N,const RealType a,const RealType b,
            VectorType & a_out, VectorType & b_out)
        {
            //TODO : make use the eigen assert facilities
            assert(a > -1);
            assert(b > -1);
            assert(a_out.rows()==b_out.rows());
            assert(a_out.rows() > 0);
            assert(N<=a_out.rows());

            a_out(0) = (b-a)/(a+b+2.);
            b_out(0) = pow(2.,(a+b+1.))*Gamma(a+1.)*Gamma(b+1.)/Gamma(a+b+2.);

            for(IndexType n=1;n<N;++n)
            {
                RealType nab = 2.*n+a+b;
                a_out(n) = (b*b - a*a)/(nab*(nab+2.));
                b_out(n) =  4.*(n+a)*(n+b)*n*(n+a+b)/(nab*nab*(nab+1.)*(nab-1.));
            }
        }

        static void r_jacobi_01(const IndexType N,const RealType a,const RealType b,
            VectorType & a_out, VectorType & b_out)
        {
            //TODO : make use the eigen assert facilities
            assert(a > -1);
            assert(b > -1);
            assert(a_out.rows()==b_out.rows());
            assert(a_out.rows() > 0);
            assert(N<=a_out.rows());

            r_jacobi(N,a, b, a_out, b_out);

            for(IndexType n=0;n<N;++n)
            {
                a_out(n)  = (1+a_out(n))/2.;
            }

            b_out(0) = b_out(0)/pow(2,a+b+1.);

            for(IndexType n=1;n<N;++n)
            {
                b_out(n) = b_out(n)/4.;
            }

        }


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
                IndexType k;
                IndexType j;
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

            //TODO : CHECK NEEDED LIKE THE ONE ON LINE 21 IN KONRAD.M
            // Do we have an approximately equal function in Eigen?
            assert( fabs(b0.sum() - (RealType) (2*N+1)) > 1e-5 );

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

        static void mpkonrad(const IndexType N,VectorType & x, VectorType & w)
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
                x(i) = 2.*x(i) - 1.;
                w(i) = 2.*w(i);
            }
        }

        static void mpkonrad15(Eigen::Array<RealType, 8, 1> & abscissaeGaussKronrod15,
            Eigen::Array<RealType, 8, 1> & weightsGaussKronrod15)
        {
            const IndexType N = 7;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod15 = xout;
            weightsGaussKronrod15 = wout;
        }

static Eigen::Array<RealType, 8, 1> mpkonrad15abscissae()
{
    const IndexType N = 7;
    VectorType x=VectorType::Zero(2*N+1);
    VectorType w=VectorType::Zero(2*N+1);

    mpkonrad(N,x,w);

    Eigen::Array<RealType, 8, 1> xout;

    IndexType ind=0;
    for(IndexType i=2*N;i>=N;--i)
    {
        xout(ind) = x(i);
        ++ind;
    }

    return xout;
}

static Eigen::Array<RealType, 8, 1> mpkonrad15weights()
{
    const IndexType N = 7;
    VectorType x=VectorType::Zero(2*N+1);
    VectorType w=VectorType::Zero(2*N+1);

    mpkonrad(N,x,w);

    Eigen::Array<RealType, 8, 1> wout;

    IndexType ind=0;
    for(IndexType i=2*N;i>=N;--i)
    {
        wout(ind) = w(i);
        ++ind;
    }

    return wout;
}

        static void mpkonrad21(Eigen::Array<RealType, 11, 1> & abscissaeGaussKronrod21,
            Eigen::Array<RealType, 11, 1> & weightsGaussKronrod21)
        {
            const IndexType N = 10;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod21 = xout;
            weightsGaussKronrod21 = wout;
        }

        static void mpkonrad31(Eigen::Array<RealType, 16, 1> & abscissaeGaussKronrod31,
            Eigen::Array<RealType, 16, 1> & weightsGaussKronrod31)
        {
            const IndexType N = 15;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod31 = xout;
            weightsGaussKronrod31 = wout;
        }


        static void mpkonrad41(Eigen::Array<RealType, 21, 1> & abscissaeGaussKronrod41,
            Eigen::Array<RealType, 21, 1> & weightsGaussKronrod41)
        {
            const IndexType N = 20;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod41 = xout;
            weightsGaussKronrod41 = wout;
        }


        static void mpkonrad51(Eigen::Array<RealType, 26, 1> & abscissaeGaussKronrod51,
            Eigen::Array<RealType, 26, 1> & weightsGaussKronrod51)
        {
            const IndexType N = 25;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod51 = xout;
            weightsGaussKronrod51 = wout;
        }


        static void mpkonrad61(Eigen::Array<RealType, 31, 1> & abscissaeGaussKronrod61,
            Eigen::Array<RealType, 31, 1> & weightsGaussKronrod61)
        {
            const IndexType N = 30;
            VectorType x=VectorType::Zero(2*N+1);
            VectorType w=VectorType::Zero(2*N+1);

            mpkonrad(N,x,w);

            Eigen::Array<RealType, N+1, 1> xout;
            Eigen::Array<RealType, N+1, 1> wout;

            IndexType ind=0;
            for(IndexType i=2*N;i>=N;--i)
            {
                xout(ind) = x(i);
                wout(ind) = w(i);
                ++ind;
            }

            abscissaeGaussKronrod61 = xout;
            weightsGaussKronrod61 = wout;
        }
    };
}

#endif //NI_KRONRADLAURIEGAUTSCHI_STB_HPP
