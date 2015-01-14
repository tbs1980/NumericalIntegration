/**
 * \file
 * \brief - This routine calculates an approximate integral to a Cauchy principal value i = integeral of f*w over (lowerLimit, upperLimit)
 *          (w(x) = 1/(x-c), (c != a, c != b), hopefully satisfying following claim for accuracy abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i))
 * \details - This function allows for the computation of a Cauchy principal value employing a Clenshaw-Curtis method automatic integerator.
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[in] c - Parameter in the weight function, c != a, c != b if c = a or c = b, the routine will end with m_errorCode = 6.
 * \param[in] desiredAbsoluteError - Absolute accuracy requested
 * \param[in] desiredRelativeError - Relative accuracy requested if  desiredAbsoluteError <= 0 and desiredRelativeError < max(50*rel.mach.acc.,0.5d-28), the routine will end with m_errorCode = 6.
 * \param[in] integral - Approximation to the integeral
 * \param[in] m_estimatedError - Estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)
 * \param[in] maxErrorIndex - pointer to the interval with largest error estimate
 * \param[in] errorMax - list(maxErrorIndex)
 * \param[in] area - sum of the integerals over the subintervals
 * \param[in] errorSum - sum of the errors over the subintervals
 * \param[in] errorBound - requested accuracy max(desiredAbsoluteError,desiredRelativeError*abs(integral))
 *
 * \returns The approximation to the integeral.
 */
 
 template <typename FunctionType>
Scalar adaptiveQuadratureForCauchyPrincipalValue(
        const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
        const Scalar desireabsoluteError = Scalar(0.), const Scalar desiredRelativeError = Scalar(0.),
        const QuadratureRule quadratureRule = 1)
{
    if ((desiredAbsoluteError <= Scalar(0.) && desiredRelativeError < Eigen::NumTraits<Scalar>::epsilon())
        || m_maxSubintervals < 1)
    {
        m_errorCode = 6;
        return Scalar(0.);
    }

    m_errorCode = 0;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    m_lowerList[0] = lowerLimit;
    m_upperList[0] = upperLimit;
    m_integralList[0] = Scalar(0.);
    list[0] = Scalar(0.);
    listIndices[0] = 0;
    listIndices[1] = 1;

    Scalar integral = 0.;
    Scalar absDiff = 0.;
    Scalar m_estimatedError = 0.;
    Scalar absIntegral = 0.;
      real a,aa,m_estimatedError,m_lowerList,area,area1,area12,area2,a1,a2,b,bb,m_upperList,b1,b2,c,r1mach,list,Eigen::NumTraits<Scalar>::epsilon(),desiredAbsoluteError,desiredRelativeError,errorBound,errorMax,error1,error2,errorSum,f,integral,m_integralList,(std::numeric_limits<Scalar>::min)()
      integer m_errorCode,listIndices,roundOff1,roundOff2,k,ruleKeye,m_numSubintervals,m_maxSubintervals,maxErrorIndex,nEval,m_numEvaluations,maxNumberOfIntegrals

    list(m_maxSubintervals);
    listIndices(m_maxSubintervals);

    m_estimatedError = 0.;

    if(c == a || c == b || (desiredAbsoluteError <= 0 && desiredRelativeError < 50. * Eigen::NumTraits<Scalar>::epsilon())
    {
        m_errorCode = 6;
        return Scalar(0.);
    }


    //first approximation to the integeral
      aa=lowerLimit;
      bb=upperLimit;

    if (lowerLimit > upperLimit)
    {
        aa=upperLimit;
        bb=lowerLimit;
    }

      ruleKeye = 1;
      call qc25c(f,alowerLimit, upperLimit,c,integral,m_estimatedError,ruleKeye,m_numEvaluations);
      
      m_numSubintervals = 1
      m_integralList(1) = integral
      list(1) = m_estimatedError
      listIndices(1) = 1
      m_lowerList(1) = a
      m_upperList(1) = b

  // Test on accuracy
      errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(integral))
      if(m_maxSubintervals == 1) m_errorCode = 1
      if(m_estimatedError < (std::min)(0.01*abs(integral),errorBound) || m_errorCode == 1) go to 70

      m_lowerList(1) = aa
      m_upperList(1) = bb
      m_integralList(1) = integral
      errorMax = m_estimatedError
      maxErrorIndex = 1
      area = integral
      errorSum = m_estimatedError
      maxNumberOfIntegrals = 1
      roundOff1 = 0
      roundOff2 = 0

  // Main loop for the integration
      do 40 m_numSubintervals = 2,m_maxSubintervals

  // Bisect the subinterval with maxNumberOfIntegrals-th largest error estimate.
        a1 = m_lowerList(maxErrorIndex)
        b1 = 0.5*(m_lowerList(maxErrorIndex)+m_upperList(maxErrorIndex))
        b2 = m_upperList(maxErrorIndex)
        if(c <= b1 && c > a1) b1 = 0.5*(c+b2)
        if(c > b1 && c < b2) b1 = 0.5*(a1+c)
        a2 = b1
        ruleKeye = 2
        call qc25c(f,a1,b1,c,area1,error1,ruleKeye,nEval)
        m_numEvaluations = m_numEvaluations+nEval
        call qc25c(f,a2,b2,c,area2,error2,ruleKeye,nEval)
        m_numEvaluations = m_numEvaluations+nEval

  // Improve previous approximations to integeral and error and test for accuracy.
        area12 = area1+area2
        erro12 = error1+error2
        errorSum = errorSum+erro12-errorMax
        area = area+area12-m_integralList(maxErrorIndex)
        if(abs(m_integralList(maxErrorIndex)-area12) < 0.1e-04*abs(area12) && erro12 >= 0.99 *errorMax && ruleKeye == 0) roundOff1 = roundOff1+1
        if(m_numSubintervals > 10 && erro12 > errorMax && ruleKeye == 0) roundOff2 = roundOff2+1
        m_integralList(maxErrorIndex) = area1
        m_integralList(m_numSubintervals) = area2
        errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(area))
        if(errorSum <= errorBound) go to 15

    // Test for roundoff error and eventually set error flag.
    if(roundOff1 >= 6 && roundOff2 > 20) m_errorCode = 2

    // Set error flag in the case that number of interval bisections exceeds m_maxSubintervals.
    if(m_numSubintervals == m_maxSubintervals) m_errorCode = 1

    // Set error flag in the case of bad integerand behaviour at a point of the integeration range.
    if((std::max)(abs(a1),abs(b2)) <= (1.+100.*Eigen::NumTraits<Scalar>::epsilon())*(abs(a2)+0.1e+04*(std::numeric_limits<Scalar>::min)()))
    {
        m_errorCode = 3;
    }

    // append the newly-created intervals to the list.
   15   if(error2 > error1) go to 20
        m_lowerList(m_numSubintervals) = a2
        m_upperList(maxErrorIndex) = b1
        m_upperList(m_numSubintervals) = b2
        list(maxErrorIndex) = error1
        list(m_numSubintervals) = error2
        go to 30
   20   m_lowerList(maxErrorIndex) = a2
        m_lowerList(m_numSubintervals) = a1
        m_upperList(m_numSubintervals) = b1
        m_integralList(maxErrorIndex) = area2
        m_integralList(m_numSubintervals) = area1
        list(maxErrorIndex) = error2
        list(m_numSubintervals) = error1

    // call subroutine qpsrt to maintain the descending ordering in the list of error estimates and select the
    // subinterval with maxNumberOfIntegrals-th largest error estimate (to be bisected next).
   30    call qpsrt(m_maxSubintervals,m_numSubintervals,maxErrorIndex,errorMax,list,listIndices,maxNumberOfIntegrals)
    // jump out of do-loop
        if(m_errorCode != 0 || errorSum <= errorBound) go to 50
   40 continue

    // compute final integral.
   50 integral = 0.
      do 60 k=1,m_numSubintervals
        integral = integral+m_integralList(k)
   60 continue
      m_estimatedError = errorSum
   70 if (aa == b) integral=-integral
  999 return
      end
