/**
 * \file
 * \brief - This routine computes modified chebsyshev moments. the k-th modified Chebyshev moment is defined as the integeral over
 *          (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev polynomial of degree k.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] alpha - parameter in the weight function w(x), alpha > (-1)
 * \param[] beta - parameter in the weight function w(x), beta > (-1)
 * \param[] ri - vector of dimension 25.  ri(k) is the integeral over (-1,1) of pow((1+x), alpha) * t(k-1,x) from k = 1 to 25.
 * \param[] rj - vector of dimension 25. rj(k) is the integeral over (-1,1) of pow((1-x), beta) * t(k-1,x) from k = 1 to 25.
 * \param[] rg - vector of dimension 25. rg(k) is the integeral over (-1,1) of pow((1+x),alpha) * log((1+x)/2) * t(k-1,x) from k = 1 to 25.
 * \param[] rh - vector of dimension 25. rh(k) is the integeral over (-1,1) of pow((1-x),beta) * log((1-x)/2) * t(k-1,x) from k = 1 to 25.
 * \param[] momentRule - input parameter indicating the modified moments to be computed
 *                       momentRule = 1 compute ri, rj
 *                                  = 2 compute ri, rj, rg
 *                                  = 3 compute ri, rj, rh
 *                                  = 4 compute ri, rj, rg, rh
 *
 * \returns The approximation to the integeral.
 */
 
      Scalar aN;
      Scalar alpha;
      Scalar beta;
      Scalar rg;
      Scalar rh;
      Scalar ri;
      Scalar rj;

      int i;
      int iMinus1;
      int momentRule;

      dimension rg(25),rh(25),ri(25),rj(25)

    // compute ri, rj, rg, rh using forward recurrence relations.
    ri(1) = pow(2.,alpha+1) / (alpha+1);
    rj(1) = pow(2.,beta+1) / (beta+1);
    ri(2) = ri(1)*alpha / (alpha+2);
    rj(2) = rj(1)*beta / (beta+2);
    aN = 2.;
    
    using std::pow;
    
    for(size_t i=2; i<25; ++i)
    {
        ri(i) = -(pow(2.,alpha+1) + aN * (aN-alpha+2) * ri(i-1)) / ((aN-1) * (aN+alpha+1));
        rj(i) = -(pow(2.,beta+1) + aN * (aN-beta+2) * rj(i-1))/((aN-1) * (aN+beta+1));
        ++aN;
    }

    if(momentRule == 1)
    {
        for (size_t i=1; i<25; i+2)
        {
            rj(i) = -rj(i)
            return;
        }
    }

    if(momentRule !=3)
    {
        rg(1) = -ri(1) / (alpha+1);
        rg(2) = -(pow(2., alpha+1) + pow(2.,alpha+1)) / (alpha+2*alpha+2)-rg(1);
        an = 2.
        aN-1 = 1.
        iMinus1 = 2

        for (size_t i=2; i<25; ++i)
        {
            rg(i) = -(an * (aN-alpha+2) * rg(iMinus1) - aN*ri(iMinus1)+(aN-1) * ri(i)) / ((aN-1) * (aN+alpha+1));
            aN++
            iMinus1 = i
        }
        if(momentRule == 2)
        {
            for (size_t i=1; i<25; i+2)
            {
                rj(i) = -rj(i)
                return;
            }
        }
    }

    rh(1) = -rj(1) / (beta+1)
    rh(2) = -(pow(2.,beta+1) + pow(2.,beta+1)) / (beta+2 * beta+2) - rh(1)
    aN = 2.;

    iMinus1 = 2
    for (size_t i=2; i<25; ++i)
    {
        rh(i) = -(an*(an - beta + 2) * rh(iMinus1) - aN * rj(iMinus1) + aN-1 * rj(i)) / (aN-1 * (aN + beta + 1));
        aN++;
        iMinus1 = i;
        rh(i) = -rh(i);
        rj(i) = -rj(i);
    }
    return;
}
