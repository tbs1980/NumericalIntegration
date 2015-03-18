#ifndef EIGEN_MONEGATO_H
#define EIGEN_MONEGATO_H

namespace Eigen{

    /**
     * \ingroup Quadrature_Module
     *
     * \class Monegato
     *
     * \brief A class for computing the Gauss/Kronrod weights.
     *
     * \tparam _Scalar floating point type
     *
     * This class is based on "Some remarks on the construction of
     * extended Gaussian quadrature rules", Giovanni Monegato,
     * Math. Comp., Vol. 32 (1978) pp. 247-252. http://www.jstor.org/stable/2006272 .
     *
     * The code is based on quadpackcpp library written by
     * Jerry Gagelman <jerry@os-scientific.org> and is distributed under
     * GNU General Public License
     *
     */
    template<typename _Scalar>
    class Monegato
    {
    public:
        typedef _Scalar Scalar;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> ScalarArrayType;
        typedef typename ScalarArrayType::Index Index;

        /**
         * \brief Compute the absolute value of a scalar
         * \param[in] x inpute scalar value
         * \return absolute value of scalar input
         */
        static Scalar abs(Scalar x)
        {
            return (x < Scalar(0)) ? -(x) : x;
        }

        /**
         * \brief Compute the Legendre polynomials and their error
         * \param[in] n degree
         * \param[in] x value of the variable
         * \param[out] err error to be returned
         * \return value of the Legendre polynomial
         *
         * This function is based on the routine gsl_sf_legendre_Pl_e
         * distributed with GSL
         */
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

        /**
         * \brief Three-term recursion identity for the Legendre derivatives
         * \param[in] n degree
         * \param[in] x value of the variable
         * \return value of the Legendre derivative
         */
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

        /**
         * \brief Compute derivatives of Chebyshev polynomials
         * \param[in] x variable
         * \param[in] n_ degree
         * \param[in] coefs Chebyshev coefficients
         * \return Derivative of Chebyshev polynomials
         */
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

        /**
         * \breif Evaluation of the Chebyshev polynomial using Clenshaw recursion
         * \param[in] x variable
         * \param[in] n_ degree
         * \param[in] coefs Chebyshev coefficients
         * \param[out] err error to be returned
         * \return value of Chebyshev polynomial
         */
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

        /**
         * \brief Computes the zeros of the Legendre polynomial
         * \param[in] m_ degree
         * \param[out] zeros the zeros of the Legendre polynomial
         */
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

        /**
         * \brief Computes coefficients of Chebyshev polynomial
         * \param[in] m_ degree
         * \param[out] coefs coefficients of Chebyshev polynomial
         */
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

        /**
         * \brief Compute the Gauss/Kronrod abscissae
         * \param[in] n_ size of Gauss-Kronrod arrays
         * \param[in] Gauss-Legendre degree
         * \param[in] coefs Chebyshev coefficients
         * \param[out] xgk_ Gauss/Kronrod abscissae
         */
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

        /**
         * \brief Compute Gauss/Kronrod weights
         * \param[in] n_ size of Gauss-Kronrod arrays
         * \param[in] m_ Gauss-Legendre degree
         * \param[in] xgk_ Gauss/Kronrod abscissae
         * \param[out] wg_ Gauss weights
         * \param[out] wgk_ Kronrod weights
         */
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

        /**
         * \brief A method ofr computing Gauss/Kronrod abscissae and weights
         * \param[in] m_ Gauss-Legendre degree
         * \param[out] xgk_ Gauss/Kronrod abscissae
         * \param[out] wgk_ Gauss/Kronrod weights
         * \param[out] xk_ Gauss abscissae
         * \param[out] wg_ Gauss weights
         *
         * Note that Gauss abscissae is not calculated here.
         */
        static void computeAbscissaeAndWeights(unsigned int m_,
            ScalarArrayType& xgk_,ScalarArrayType& wgk_,
            ScalarArrayType& xk_,ScalarArrayType& wg_)
        {
            //const Index m_ = 2*nNodes;
            const unsigned int n_ = m_ + 1;

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

}//namespace Eigen

#endif //EIGEN_MONEGATO_H
