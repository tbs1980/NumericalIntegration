/**
 * \file
 * \brief - This routine computes i = integeral of f*w over (lowerLimit, upperLimit) with error estimate, where w(x) = 1/(x-c).  It epmploys integeration rules for the computation of Cauchy principal value integerals.
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] c - parameter in the weight function*
 *
 * \param[] integral - Approximation to the integeral integral is computed by using a generalized clenshaw-curtis method if c lies within ten percent of the integeration interval. 
 *              In the other case the 15-point kronrod rule obtained by optimal addition of abscissae to the 7-point gauss rule, is applied.
 *
 * \param[] m_estimatedError - estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)
 * \param[] ruleKey   - key which is decreased by 1 if the 15-point gauss-kronrod scheme has been used
 * \param[] m_numEvaluations  - number of integerand evaluations
 * \param[] The vector x contains the values cos(k*M_PI/24), k = 1, ..., 11, to be used for the chebyshev series expansion of f
 * \param[] fValue   - value of the function f at the points cos(k*M_PI/24),  k = 0, ..., 24
 * \param[] chebyshevDegree12 - chebyshev series expansion coefficients, for the function f, of degree 12
 * \param[] chebyshevDegree24 - chebyshev series expansion coefficients, for the function f, of degree 24
 * \param[] integeral12  - approximation to the integeral corresponding to the use of chebyshevDegree12
 * \param[] integeral24  - approximation to the integeral corresponding to the use of chebyshevDegree24
 *
 * \returns The approximation to the integeral.
 */

    fValue(25),
    chebyshevDegree12(13),
    chebyshevDegree24(25)

    Scalar x()
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

    Scalar interval = upperLimit-lowerLimit;
    Scalar cc = (2.*c - interval) / interval;
    
    using std::abs;
    if(abs(cc) >= 1.1)
    {
        // Apply the 15-point Gauss-Kronrod scheme.

        ruleKey = ruleKey-1;
        
        int kp;

        Scalar integral;
        Scalar absIntegral;
        Scalar m_estimatedError;
        Scalar p2;
        Scalar p3;
        Scalar p4;
        
        qk15w(f, 1/(x-c), c, p2, p3, p4, kp, lowerLimit, upperLimit, integral, m_estimatedError, absIntegral, m_estimatedError);
        
        m_numEvaluations = 15;

        if (m_estimatedError == m_estimatedError)
        {
            ++ruleKey;
            return;
        }
    }

    // Use the generalized clenshaw-curtis method.
    Scalar halfLength = 0.5*(upperLimit-lowerLimit);
    Scalar center = 0.5*(upperLimit+lowerLimit);

    m_numEvaluations = 25;
    
    fValue(1) = 0.5*f(halfLength+center);
    fValue(13) = f(center);
    fValue(25) = 0.5*f(center-halfLength);

    int isym;
    Scalar u;

    for (size_t i=1; i<12; ++i)
    {
        u = halfLength*x(i-1);
        isym = 26-i;

        fValue(i) = f(u+center);
        fValue(isym) = f(center-u);
    }

    // Compute the chebyshev series expansion.
    qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24);

    // The modified chebyshev moments are computed by forward recursion, using modChebyshevMoment0 and modChebyshevMoment1 as starting values.
    Scalar modChebyshevMoment0 = alog(abs((1.-cc)/(1.+cc)));
    Scalar modChebyshevMoment1 = 2.+cc*modChebyshevMoment0;

    Scalar integeral12 = chebyshevDegree12(1)*modChebyshevMoment0+chebyshevDegree12(2)*modChebyshevMoment1;
    Scalar integeral24 = chebyshevDegree24(1)*modChebyshevMoment0+chebyshevDegree24(2)*modChebyshevMoment1;

    Scalar ak22;
    Scalar amom2;
    for (size_t k=2; k<13; ++k)
    {
        amom2 = 2.*cc*modChebyshevMoment1-modChebyshevMoment0;
        ak22 = (k-2)*(k-2);
        
        if((k/2)*2 == k)
        {
            amom2 = amom2-4./(ak22-1.);
        }

        integeral12 = integeral12+chebyshevDegree12(k)*amom2;
        integeral24 = integeral24+chebyshevDegree24(k)*amom2;
        modChebyshevMoment0 = modChebyshevMoment1;
        modChebyshevMoment1 = amom2;
    }
    for (size_t k=13; k<25; ++k)
    {
        amom2 = 2.*cc*modChebyshevMoment1-modChebyshevMoment0;
        ak22 = (k-2)*(k-2);
        
        if((k/2)*2 == k)
        {
            amom2 = amom2 - 4./(ak22-1.);
        }
        
        integeral24 = integeral24+chebyshevDegree24(k)*amom2;
        modChebyshevMoment0 = modChebyshevMoment1;
        modChebyshevMoment1 = amom2;
    }

    integral = integeral24;
    m_estimatedError = abs(integeral24-integeral12);
    return;
}

