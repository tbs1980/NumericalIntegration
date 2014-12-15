/**
* subroutine qk15w(f,w,p1,p2,p3,p4,kp,lowerLimit, upperLimit,integral,m_estimatedError,absIntegral,m_estimatedError)
* purpose - To compute i = integeral of f*w over (lowerLimit, upperLimit), with error estimate
*           j = integeral of abs(f*w) over (lowerLimit, upperLimit)
* description - integeration rules
*
*   p1, p2, p3, p4 - Parameters in the weight function*
*
*   kp - Key for indicating the type of weight function
*
*   absIntegral - Approximation to the integeral of abs(f)
*
*   m_estimatedError - Approximation to the integeral of abs(f-i/(upperLimit-lowerLimit))
*
*   The abscissae and weights are given for the interval (-1,1).
*   Because of symmetry only the positive abscissae and their corresponding weights are given.
*
*   xgk - Abscissae of the 15-point gauss-kronrod rule xgk(2), xgk(4), ... abscissae of the 7-point Gauss rule
*           xgk(1), xgk(3), ... abscissae which are optimally added to the 7-point gauss rule
*
*   wgk - Weights of the 15-point gauss-kronrod rule
*
*   wg - Weights of the 7-point gauss rule
*
*   center  - mid point of the interval
*
*   halfLength  - half-length of the interval
*
*   absc*  - abscissa
*
*   fValue*  - function value
*
*   resultGauss   - integral of the 7-point gauss formula
*
*   resultKronrod   - integral of the 15-point kronrod formula
*
*   resultKronrodh  - Approximation to the mean value of f*w over (lowerLimit, upperLimit),
*            i.e. to i/(upperLimit-lowerLimit)
*/
                       
    Scalar a;
    Scalar absc;
    Scalar absc1;
    Scalar absc2;
    Scalar m_estimatedError;
    Scalar b;
    Scalar center;
    Scalar dhalfLength;
    Scalar r1mach;
    Scalar f;
    Scalar fc;
    Scalar fsum;
    Scalar fValue1;
    Scalar fValue2;
    Scalar fv1;
    Scalar fv2;
    Scalar halfLength;
    Scalar p1;
    Scalar p2;
    Scalar p3;
    Scalar p4;
    Scalar absIntegral;
    Scalar m_estimatedError;
    Scalar resultGauss;
    Scalar resultKronrod;
    Scalar resultKronrodh;
    Scalar integral;
    Scalar w;
    Scalar wg;
    Scalar wgk;
    Scalar xgk;

    int j;
    int jtw;
    int jtwm1;
    int kp;

    fv1(7);
    fv2(7);
    xgk(8);
    wgk(8);
    wg(4)

    xgk
    {
        0.9914553711208126,
        0.9491079123427585,
        0.8648644233597691,
        0.7415311855993944,
        0.5860872354676911,
        0.4058451513773972,
        0.2077849550078985,
        0.0000000000000000
    }.finished();

    wgk
    {
        0.02293532201052922,
        0.06309209262997855,
        0.1047900103222502,
        0.1406532597155259,
        0.1690047266392679,
        0.1903505780647854,
        0.2044329400752989,
        0.2094821410847278
    }.finished();

    wg
    {
        0.1294849661688697 ,
        0.2797053914892767 ,
        0.3818300505051889 ,
        0.4179591836734694
    }.finished();

    using std::abs;
    center = 0.5*(lowerLimit+upperLimit)
    halfLength = 0.5*(upperLimit-lowerLimit)
    dhalfLength = abs(halfLength)

    // Compute the 15-point kronrod approximation to the integeral, and estimate the error.
      fc = f(center)*w(center,p1,p2,p3,p4,kp)
      resultGauss = wg(4)*fc
      resultKronrod = wgk(8)*fc
      absIntegral = abs(resultKronrod)
      do 10 j=1,3
        jtw = j*2
        absc = halfLength*xgk(jtw)
        absc1 = center-absc
        absc2 = center+absc
        fValue1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fValue2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtw) = fValue1
        fv2(jtw) = fValue2
        fsum = fValue1+fValue2
        resultGauss = resultGauss+wg(j)*fsum
        resultKronrod = resultKronrod+wgk(jtw)*fsum
        absIntegral = absIntegral+wgk(jtw)*(abs(fValue1)+abs(fValue2))
   10 continue
      do 15 j=1,4
        jtwm1 = j*2-1
        absc = halfLength*xgk(jtwm1)
        absc1 = center-absc
        absc2 = center+absc
        fValue1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
        fValue2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
        fv1(jtwm1) = fValue1
        fv2(jtwm1) = fValue2
        fsum = fValue1+fValue2
        resultKronrod = resultKronrod+wgk(jtwm1)*fsum
        absIntegral = absIntegral+wgk(jtwm1)*(abs(fValue1)+abs(fValue2))
   15 continue
      resultKronrodh = resultKronrod*0.5
      m_estimatedError = wgk(8)*abs(fc-resultKronrodh)
      do 20 j=1,7
        m_estimatedError = m_estimatedError+wgk(j)*(abs(fv1(j)-resultKronrodh)+abs(fv2(j)-resultKronrodh))
   20 continue
      integral = resultKronrod*halfLength
      absIntegral = absIntegral*dhalfLength
      m_estimatedError = m_estimatedError*dhalfLength
      m_estimatedError = abs((resultKronrod-resultGauss)*halfLength)
      if(m_estimatedError != 0. && m_estimatedError != 0.) m_estimatedError = m_estimatedError*(std::min)(1.,(0.2e+03*m_estimatedError/m_estimatedError)**1.5 )
      if(absIntegral > (std::numeric_limits<Scalar>::min)()/(5.*Eigen::NumTraits<Scalar>::epsilon())) m_estimatedError = (std::max)((Eigen::NumTraits<Scalar>::epsilon()*5.)*absIntegral,m_estimatedError)
      return
      end
