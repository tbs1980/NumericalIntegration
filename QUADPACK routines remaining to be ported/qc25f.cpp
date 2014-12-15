/**
* subroutine qc25f(f,lowerLimit, upperLimit,omega,integer,nrMoments,maxp1,keyMomentsComputed,integral,m_estimatedError,m_numEvaluations,absIntegral,m_estimatedError,momentsComputed,chebyshevMoments)
* keywords  Integeration rules for functions with cos or sin factor, clenshaw-curtis, gauss-kronrod
* purpose - Computes the integeral i=integeral of f(x) over (lowerLimit, upperLimit) where 
*           w(x) = cos(omega*x) or (wx)=sin(omega*x) and to compute j=integeral of abs(f) 
*           over (lowerLimit, upperLimit). for smallestIntervalLength value of omega or smallestIntervalLength intervals 
*           (lowerLimit, upperLimit) 15-point Gauss-Kronrod rule used. 
*           Otherwise generalized Clenshaw-Curtis is used
*
* description - integeration rules for functions with cos or sin factor
*
*   omega  - A parameter in the weight function
*
*   integer - indicates which weight function is to be used: if integer = 1, w(x) = cos(omega*x), and  if integer = 2, w(x) = sin(omega*x)
*
*   nrMoments  - The length of interval (lowerLimit, upperLimit) is equal to the length of the original integeration interval divided by
*                2^nrMoments (we suppose that the routine is used in an adaptive integeration process, otherwise set nrMoments = 0).
*                nrMoments must be zero at the first call.
*
*   maxp1  - Gives an upper finiteBound on the number of chebyshev moments which can be stored, 
*            i.e. for the intervals of lengths abs(bupperLimit-lowerLimita)*2**(-l), for l = 0 to maxp1-2.
*
*   keyMomentsComputed - key which is one when the moments for the current interval have been computed
*
*   momentsComputed - for each interval length we need to compute the Chebyshev moments. momentsComputed counts the number of
*                     intervals for which these moments have already been computed. if nrMoments < momentsComputed or keyMomentsComputed = 1, the
*                     Chebyshev moments for the interval (lowerLimit, upperLimit) have already been computed and stored, otherwise we
*                     compute them and we increase momentsComputed.
*
*   chebyshevMoments -  array of dimension at least (maxp1,25) containing the modified chebyshev moments for the first momentsComputed interval lengths
*
*   fValue   - value of the function f at the points (upperLimit-lowerLimit)*0.5*cos(k*M_PI/12) + (upperLimit+lowerLimit)*0.5, from k = 0 to 24
*
*   chebyshevDegree12 - coefficients of the chebyshev series expansion of degree 12, for the function f, in the interval (lowerLimit, upperLimit)
*
*   chebyshevDegree24 - coefficients of the chebyshev series expansion of degree 24, for the function f, in the interval (lowerLimit, upperLimit)
*
*   integralCosineChebyshevDegree12 - approximation to the integeral of 
*           cos(0.5*(upperLimit-lowerLimit)*omega*x)*f(0.5*(upperLimit-lowerLimit)*x+0.5*(upperLimit+lowerLimit))
*           over (-1,+1), using the chebyshev series expansion of degree 12
*
*   integralCosineChebyshevDegree24 - approximation to the same integeral, using the chebyshev series expansion of degree 24
*
*   integralSineChebyshevDegree12 - the analogue of integralCosineChebyshevDegree12 for the sine
*
*   integralSineChebyshevDegree24 - the analogue of integralCosineChebyshevDegree24 for the sine
*/

    Scalar ac;
    Scalar an;
    Scalar an2;
    Scalar as;
    Scalar asap;
    Scalar ass;
    
    Scalar center;
    Scalar chebyshevMoments;

    Scalar cosineCenterOmega;
    Scalar sineCenterOmega;

    Scalar cosineIntegrandParam;
    Scalar sineIntegrandParam;
    
    Scalar estc;
    Scalar ests;

    Scalar d;
    Scalar d1;
    Scalar d2;

    Scalar f;
    Scalar fValue;
    Scalar halfLength;
    
    Scalar omega;
    Scalar integrandParam;

    Scalar par2;
    Scalar par22;
    Scalar p2;
    Scalar p3;
    Scalar p4;

    
    Scalar integralCosineChebyshevDegree12;
    Scalar integralCosineChebyshevDegree24;
    Scalar integralSineChebyshevDegree12;
    Scalar integralSineChebyshevDegree24;

    int integer;
    int isym;

    int keyMomentsComputed;
    int maxp1;
    int momentsComputed;
    int m_numEvaluations;
    int noequ;
    int noeq1;
    int nrMoments;

    chebyshevMoments(maxp1, 25);
    chebyshevDegree12(13);
    chebyshevDegree24(25);
    d(25);
    d1(25);
    d2(25);
    fValue(25);
    v(28);
    x(11);

    // The vector x contains the values cos(k*M_PI/24) for k = 1 to 11, to be used for the chebyshev expansion of f
    Scalar x
    {
        0.9914448613738104,
        0.9659258262890683,
        0.9238795325112868,
        0.8660254037844386,
        0.7933533402912352,
        0.7071067811865475,
        0.6087614290087206,
        0.5000000000000000,
        0.3826834323650898,
        0.2588190451025208,
        0.1305261922200516
    }.finished()

      center = 0.5 * (upperLimit+lowerLimit);
      halfLength = 0.5 * (upperLimit-lowerLimit);
      integrandParam = omega * halfLength;

    if (integer == 1)
    {
        qwgtf = cos(omega*x);
    }
    else
    {
        qwgtf = sin(omega*x);
    }

    using std::abs;

    // Compute the integeral using the 15-point gauss-kronrod formula if the value of the parameter in the integerand is smallestIntervalLength.
    if(abs(integrandParam) <= 2.)
    {
        qk15w(f,qwgtf,omega,p2,p3,p4,integer,lowerLimit, upperLimit,integral,m_estimatedError,absIntegral,m_estimatedError);
        m_numEvaluations = 15;
        return integral;
    }

    // Compute the integeral using the generalized Clenshaw-Curtis method.
    using std::cos;
    cosineCenterOmega = halfLength*cos(center*omega);
    sineCenterOmega = halfLength*sin(center*omega);
    m_estimatedError = (std::numeric_limits<Scalar>::max)();
    m_numEvaluations = 25;

    // Check whether the Chebyshev moments for this interval have already been computed.
    if(nrMoments >= momentsComputed || keyMomentsComputed != 1)
    {
        // Compute a new set of chebyshev moments.
        m = momentsComputed+1;
        par2 = integrandParam*integrandParam;
        par22 = par2 + 2.;
        sineIntegrandParam = sin(integrandParam);
        cosineIntegrandParam = cos(integrandParam);

        // Compute the chebyshev moments with respect to cosine.
        v(1) = 2. * sineIntegrandParam / integrandParam;
        v(2) = (8. * cosineIntegrandParam + (par2 + par2 - 8.) * sineIntegrandParam / integrandParam)/par2;
        v(3) = (32.*(par2 - 12.)*cosineIntegrandParam + (2.*((par2-80.)*par2 + 192.) * sineIntegrandParam) / integrandParam) / (par2*par2);
        ac = 8.*cosineIntegrandParam;
        as = 24.*integrandParam*sineIntegrandParam;

        if(abs(integrandParam) <= 24.)
        {
            // Compute the chebyshev moments as the solutions of a finiteBoundary value problem with 1 initial value (v(3)) and 1 end value (computed using an asymptotic formula).
            noequ = 25;
            noeq1 = noequ-1;
            an = 6.;
            
            for ( size_t k=1; k<noeq1; ++k)
            {
                an2 = an*an;
                d(k) = -2.*(an2-4.)*(par22-an2-an2);
                d2(k) = (an-1.)*(an-2.)*par2;
                d1(k+1) = (an+3.)*(an+4.)*par2;
                v(k+3) = as-(an2-4.)*ac;
                an = an+2.;
            }

            an2 = an*an;
            d(noequ) = -2.*(an2-4.)*(par22-an2-an2);
            v(noequ+3) = as-(an2-4.)*ac;
            v(4) = v(4)-0.56e+02*par2*v(3);
            ass = integrandParam*sineIntegrandParam;
            asap = (((((210.*par2-1.)*cosineIntegrandParam-(105.*par2-63.)*ass)/an2-(1.-15.*par2)*cosineIntegrandParam+15.*ass)/an2-cosineIntegrandParam+3.*ass)/an2-cosineIntegrandParam)/an2;
            v(noequ+3) = v(noequ+3)-2.*asap*par2*(an-1.)*(an-2.);

            // Solve the tridiagonal system by means of Gaussian elimination with partial M_PIvoting.
            sgtsl(noequ,d1,d,d2,v(4),m_errorCodes);
        }
        else
        {
            // Compute the chebyshev moments by means of forward recursion.
            an = 4.;

            for (size_t i=4; i<13; ++i)
            {
                an2 = an*an;
                v(i) = ((an2-4.)*(2.*(par22-an2-an2)*v(i-1)-ac)+as-par2*(an+1.)*(an+2.)*v(i-2))/(par2*(an-1.)*(an-2.));
                an = an+2.;
            } 
        }


        for (size_t j=1; j<13; ++j)
        {
            chebyshevMoments(m,2*j-1) = v(j);
        }

        // Compute the chebyshev moments with respect to sine.
        v(1) = 2.*(sineIntegrandParam-integrandParam*cosineIntegrandParam)/par2;
        v(2) = (18. - 48./par2) * sineIntegrandParam / par2 + (-2. + 48. / par2) * cosineIntegrandParam / integrandParam;
        ac = -2.4 * integrandParam * cosineIntegrandParam;
        as = -0.8 * sineIntegrandParam;

        if(abs(integrandParam) <= 2.4)
        {
            // Compute the chebyshev moments as the solutions of a finiteBoundary value problem with 1 initial value (v(2)) and 1 end value (computed using an asymptotic formula).
            an = 5.;

            for (size_t k=1; k<noeq1; ++k)
            {
                an2 = an*an;
                d(k) = -2.*(an2-4.)*(par22-an2-an2);
                d2(k) = (an-1.)*(an-2.)*par2;
                d1(k+1) = (an+3.)*(an+4.)*par2;
                v(k+2) = ac+(an2-4.)*as;
                an = an+2.;
            }

            an2 = an*an;
            d(noequ) = -2.*(an2-4.)*(par22-an2-an2);
            v(noequ+2) = ac+(an2-4.)*as;
            v(3) = v(3)-42.*par2*v(2);
            ass = integrandParam*cosineIntegrandParam;
            asap = (((((105.*par2-63.)*ass+(210.*par2-1.)*sineIntegrandParam)/an2+(15.*par2-1.)*sineIntegrandParam-15.*ass)/an2-3.*ass-sineIntegrandParam)/an2-sineIntegrandParam)/an2;
            v(noequ+2) = v(noequ+2)-2.*asap*par2*(an-1.)*(an-2.);

            // Solve the tridiagonal system by means of gaussian elimination with partial M_PIvoting.
            call sgtsl(noequ,d1,d,d2,v(3),m_errorCodes);
        }
        else
        {
            // Compute the chebyshev moments by means of forward recursion.
            an = 3.;

            for size_t i=2; i<12; ++i)
            {
                an2 = an*an;
                v(i) = ((an2 - 4.) * (2. * (par22 - an2 - an2) * v(i-1) + as) + ac - par2 * (an + 1.) * (an + 2.) * v(i - 2)) / (par2 * (an - 1.) * (an - 2.));
                an += 2.;
            }
        }
    
        for (size_t j=1;j<12; ++j);
        {
            chebyshevMoments(m, 2*j) = v(j);
        }
    }
    
    if (nrMoments < momentsComputed)
    {
        m = nrMoments + 1;
    }
    
    if (momentsComputed < maxp1-1 && nrMoments >= momentsComputed)
        {
            ++momentsComputed;
        }

    // Compute the coefficients of the chebyshev expansions of degrees 12 and 24 of the function f.
    fValue(1) = 0.5 * f(center + halfLength);
    fValue(13) = f(center);
    fValue(25) = 0.5 * f(center - halfLength);

    for (size_t i=1; i<12; ++i)
    {
        isym = 26 - i;
        fValue(i) = f(halfLength * x(i - 1) + center);
        fValue(isym) = f(center-halfLength * x(i - 1));
    }

    Scalar chebyshevDegree12;
    Scalar chebyshevDegree24;
    qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24);

    // Compute the integeral and error estimates.
    integralCosineChebyshevDegree12 = chebyshevDegree12(13) * chebyshevMoments(m, 13);
    integralSineChebyshevDegree12 = 0.;
    k = 11

    for (size_t j=0; j<6; ++j)
    {
        integralCosineChebyshevDegree12 = integralCosineChebyshevDegree12+chebyshevDegree12(k) * chebyshevMoments(m, k);
        integralSineChebyshevDegree12 = integralSineChebyshevDegree12+chebyshevDegree12(k+1) * chebyshevMoments(m, k+1);
        k -=2;
    }

    integralCosineChebyshevDegree24 = chebyshevDegree24(25) * chebyshevMoments(m, 25);
    integralSineChebyshevDegree24 = 0.;
    absIntegral = abs(chebyshevDegree24(25));
    k = 23;

    for (size_t j=0; j<12; ++j)
    {
        integralCosineChebyshevDegree24 = integralCosineChebyshevDegree24 + chebyshevDegree24(k) * chebyshevMoments(m, k);
        integralSineChebyshevDegree24 = integralSineChebyshevDegree24 + chebyshevDegree24(k+1) * chebyshevMoments(m, k+1);
        absIntegral = abs(chebyshevDegree24(k)) + abs(chebyshevDegree24(k+1));
    }

    estc = abs(integralCosineChebyshevDegree24 - integralCosineChebyshevDegree12);
    ests = abs(integralSineChebyshevDegree24 - integralSineChebyshevDegree12);
    absIntegral = absIntegral * abs(halfLength);

    if(integer == 2)
    {
        integral = cosineCenterOmega * integralSineChebyshevDegree24 + sineCenterOmega * integralCosineChebyshevDegree24;
        m_estimatedError = abs(cosineCenterOmega * ests) + abs(sineCenterOmega * estc);
    }
    else
    {
        integral = cosineCenterOmega * integralCosineChebyshevDegree24 - sineCenterOmega * integralSineChebyshevDegree24;
        m_estimatedError = abs(cosineCenterOmega * estc) + abs(sineCenterOmega * ests);
    }
    
    return integral;
}
