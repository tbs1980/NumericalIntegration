/**
 * \file
 * \brief - This routine calculates an approximation integral to a given definite integeral i = integeral of 
 *           f over (lowerLimit, upperLimit),hopefully satisfying following claim for accuracy 
 *           abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i)).
 *           Break points of the integeration interval, where local difficulties of the integerand may
 *           occur(e.g. singularities,discontinuities),provided by user.
 *
 * \details - This routine is a general-purpose, globally adaptive, extrapolative, automatic integerator for functions with singularities occuring at known (user specified) points over a definite integeral.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[in] points - Vector of dimension singularityPointsPlus2, the first (singularityPointsPlus2-2) elements of which are the user provided break points.  If these points do not constitute an ascending sequence there will be an automatic sorting.
 * \param[in] desiredAbsoluteError - Absolute accuracy requested
 * \param[in] desiredRelativeError - Relative accuracy requested. If  desiredAbsoluteError <= 0 and desiredRelativeError < max(50*rel.mach.acc.,0.5d-28), the routine will end with m_errorCode = 6.
 * \param[in] integral - Approximation to the integeral
 * \param[in] m_estimatedError - Estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)
 * \param[in] singularityPoints - Vector of dimension at least singularityPointsPlus2, containing the integeration m_maxSubintervalss and the break points of the interval in ascending sequence.
 * \param[in] singularityPointsPlus2 - Number equal to two more than the number of user-supplied break points within the integeration range, singularityPointsPlus2 >= 2.  If singularityPointsPlus2 < 2, the routine will end with m_errorCode = 6.
 * \param[in] level - Vector of dimension at least m_maxSubintervals, containing the subdivision levels of the subinterval, i.e. if (alowerLimit, upperLimitb) is a subinterval of (p1,p2) where p1 as well as p2 is a user-provided break point or integeration m_maxSubintervals, then (alowerLimit, upperLimitb) has level l if abs(bupperLimit-lowerLimita) = abs(p2-p1)*2**(-l).
 * \param[in] ndin - Vector of dimension at least singularityPointsPlus2, after first integeration over the intervals (singularityPoints(i)),singularityPoints(i+1), from i = 0 to singularityPointsPlus2-2, the error estimates over some of the intervals may have been increased artificially, in order to put their subdivision forward. if this happens for the subinterval numbered k, ndin(k) is put to 1, otherwise ndin(k) = 0.
 * \param[in] maxErrorIndex - pointer to the interval with largest error estimate
 * \param[in] errorMax - list(maxErrorIndex)
 * \param[in] errorLast - error on the interval currently subdivided (before that subdivision has taken place)
 * \param[in] area      - sum of the integerals over the subintervals
 * \param[in] errorSum    - sum of the errors over the subintervals
 * \param[in] errorBound - requested accuracy max(desiredAbsoluteError,desiredRelativeError*abs(integral))
 * \param[in] numberOfExtrapolations - number of calls to the extrapolationPerformed routine
 * \param[in] numberElementsList2 - Number of elements in m_integralList2.  If an appropriate approximation to the compounded integeral has been obtained, it is put in m_integralList2(numberElementsList2) after numberElementsList2 has been increased by one.
 * \param[in] errorLargest - sum of the errors over the intervals larger than the smallestIntervalLengthest interval considered up to now
 * \param[in] extrapolationPerformed - logical variable denoting that the routine is attempting to perform extrapolationPerformed. i.e. before subdividing the smallestIntervalLengthest interval we try to decrease the value of errorLargest.
 * \param[in] extrapolationAllowed - logical variable denoting that extrapolationPerformed is still allowed (true-value)
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForSingularitiesAtKnownPoints(
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


    Scalar absoluteEpsilon;
    Scalar area;
    Scalar area1;
    Scalar area12;
    Scalar area2;
    Scalar errorCorrection;
    Scalar absDiff;
    Scalar absDiff1;
    Scalar absDiff2;
    Scalar dResult;
    Scalar r1mach;
    Scalar list;
    Scalar desiredAbsoluteError;
    Scalar desiredRelativeError;
    Scalar errorLargest;
    Scalar errorLast;
    Scalar errorBound;
    Scalar errorMax;
    Scalar error1;
    Scalar erro12;
    Scalar error2;
    Scalar errorSum;
    Scalar errorTest;
    Scalar f;
    Scalar points;
    Scalar singularityPoints;
    Scalar approximateIntegralResult;
    Scalar absIntegral;
    Scalar resultEpsilon;
    Scalar integral;
    Scalar res3la;
    Scalar sign;
    Scalar temp;

    
    int id;
    int ind1;
    int ind2;
    int ip1;
    int roundOff1;
    int roundOff2;
    int roundOff3;
    int lowerBound;
    int upperBound;
    int keySign;
    int ktmin;
    int currentLevel;
    int level;
    int levmax;
    int maxErrorIndex;
    int ndin;
    int nint;
    int nintp1;
    int nsingularityPoints;
    int singularityPointsPlus2;
    int numberOfExtrapolations;
    int maxNumberOfIntegrals;
    int numberElementsList2
    
    dimension level(m_maxSubintervals);
    ndin(singularityPointsPlus2);
    points(singularityPointsPlus2);
    singularityPoints(singularityPointsPlus2);
    res3la(3);

    m_errorCode = 0;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    integral = 0.;
    m_estimatedError = 0.;
    m_lowerList(1) = lowerLimit;
    m_upperList(1) = upperLimit;
    m_integralList(1) = 0.;
    list(1) = 0.;
    listIndices(1) = 0;
    level(1) = 0;
    nsingularityPoints = singularityPointsPlus2-2;

    if(singularityPointsPlus2 < 2 || m_maxSubintervals <= nsingularityPoints || (desiredAbsoluteError <= 0. && desiredRelativeError < (std::max)(5.*Eigen::NumTraits<Scalar>::epsilon(),0.5e-14)))
    {
        m_errorCode = 6;
        return Scalar(0.);
    }

    //if any break points are provided, sort them into an ascending sequence.
    sign = 1.0;
    if(lowerLimit > upperLimit)
    {
        sign = -1.0;
    }

    singularityPoints(1) = (std::min)(lowerLimit, upperLimit)
    if(nsingularityPoints != 0)
    {
        for (size_t i=0; i<nsingularityPoints; ++i)
        {
            singularityPoints(i+1) = points(i)
        }
    }

   15 singularityPoints(nsingularityPoints+2) = (std::max)(lowerLimit, upperLimit);
      nint = nsingularityPoints+1
      a1 = singularityPoints(1)
      if(nsingularityPoints == 0) go to 40
      nintp1 = nint+1
      do 20 i = 1,nint
        ip1 = i+1
        do 20 j = ip1,nintp1
          if(singularityPoints(i) <= singularityPoints(j)) go to 20
          temp = singularityPoints(i)
          singularityPoints(i) = singularityPoints(j)
          singularityPoints(j) = temp
   20 continue
      if(singularityPoints(1) != (std::min)(lowerLimit, upperLimit) || singularityPoints(nintp1) != (std::max)(lowerLimit, upperLimit)) m_errorCode = 6
      if(m_errorCode == 6) go to 999

    // Compute first integeral and error approximations.
   40 absIntegral = 0.
      do 50 i = 1,nint
        b1 = singularityPoints(i+1)
        call qk21(f,a1,b1,area1,error1,absDiff,approximateIntegralResult)
        m_estimatedError = m_estimatedError+error1
        integral = integral+area1
        ndin(i) = 0
        if(error1 == approximateIntegralResult && error1 != 0.) ndin(i) = 1
        absIntegral = absIntegral+absDiff
        level(i) = 0
        list(i) = error1
        m_lowerList(i) = lowerLimit1
        m_upperList(i) = upperLimit1
        m_integralList(i) = area1
        listIndices(i) = i
        lowerLimit1 = upperLimit1
   50 continue
      errorSum = 0.
      do 55 i = 1,nint
        if(ndin(i) == 1) list(i) = m_estimatedError
        errorSum = errorSum+list(i)
   55 continue

    // Test on accuracy.
      m_numSubintervals = nint
      m_numEvaluations = 21*nint
      dResult = abs(integral)
      errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*dResult)
      if(m_estimatedError <= 100.*Eigen::NumTraits<Scalar>::epsilon()*absIntegral && m_estimatedError > errorBound) m_errorCode = 2
      if(nint == 1) go to 80
      do 70 i = 1,nsingularityPoints
        lowerBound = i+1
        ind1 = listIndices(i)
        do 60 j = lowerBound,nint
          ind2 = listIndices(j)
          if(list(ind1) > list(ind2)) go to 60
          ind1 = ind2
          k = j
   60   continue
        if(ind1 == listIndices(i)) go to 70
        listIndices(k) = listIndices(i)
        listIndices(i) = ind1
   70 continue
      if(m_maxSubintervals < singularityPointsPlus2) m_errorCode = 1
   80 if(m_errorCode != 0 || m_estimatedError <= errorBound) go to 999

      m_integralList2(1) = integral
      maxErrorIndex = listIndices(1)
      errorMax = list(maxErrorIndex)
      area = integral
      maxNumberOfIntegrals = 1
      numberOfExtrapolations = 0
      numberElementsList2 = 1
      ktmin = 0;
      bool extrapolationPerformed = false;
      bool extrapolationAllowed = true;
      errorLargest = errorSum
      errorTest = errorBound
      levmax = 1
      roundOff1 = 0;
      roundOff2 = 0;
      roundOff3 = 0;
      m_errorCode = 0;
      m_estimatedError = (std::numeric_limits<Scalar>::max)();
      keySign = -1
      if(dResult >= (1.-5.*Eigen::NumTraits<Scalar>::epsilon())*absIntegral) keySign = 1

    // Main loop for the integration
    do 160 m_numSubintervals = singularityPointsPlus2,m_maxSubintervals

    // Bisect the subinterval with the maxNumberOfIntegrals-th largest error estimate.
        currentLevel = level(maxErrorIndex)+1
        lowerLimit1 = m_lowerList(maxErrorIndex)
        upperLimit1 = 0.5*(m_lowerList(maxErrorIndex)+m_upperList(maxErrorIndex))
        lowerLimit2 = upperLimit1
        upperLimit2 = m_upperList(maxErrorIndex)
        errorLast = errorMax
        call qk21(f,a1,b1,area1,error1,approximateIntegralResult,absDiff1)
        call qk21(f,a2,b2,area2,error2,approximateIntegralResult,absDiff2)

    // Improve previous approximations to integeral and error and test for accuracy.
        m_numEvaluations = m_numEvaluations+42
        area12 = area1+area2
        erro12 = error1+error2
        errorSum = errorSum+erro12-errorMax
        area = area+area12-m_integralList(maxErrorIndex)
        if(absDiff1 == error1 || absDiff2 == error2) go to 95
        if(abs(m_integralList(maxErrorIndex)-area12) > 0.1e-04*abs(area12) || erro12 < 0.99 *errorMax) go to 90
        if(extrapolationPerformed) roundOff2 = roundOff2+1
        if(.not.extrapolationPerformed) roundOff1 = roundOff1+1
   90   if(m_numSubintervals > 10 && erro12 > errorMax) roundOff3 = roundOff3+1
   95   level(maxErrorIndex) = currentLevel
        level(m_numSubintervals) = currentLevel
        m_integralList(maxErrorIndex) = area1
        m_integralList(m_numSubintervals) = area2
        errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(area))


    // Set error flag in the case that the number of subintervals equals m_maxSubintervals.
    if(m_numSubintervals == m_maxSubintervals)
    {
        m_errorCode = 1;
    }

    // Test for roundoff error and eventually set error flag.
    if(roundOff1+roundOff2 >= 10 || roundOff3 >= 20)
    {
        m_errorCode = 2;
    }

    if(roundOff2 >= 5)
    {
        m_errorCode = 3;
    }


    // Set error flag in the case of bad integerand behaviour at a point of the integeration range
    if((std::max)(abs(a1),abs(b2)) <= (1.+100.*Eigen::NumTraits<Scalar>::epsilon())*(abs(a2)+0.1e+04*(std::numeric_limits<Scalar>::min)()))
    {
        m_errorCode = 4;
    }
        // Append the newly-created intervals to the list.
        if (error2 > error1)
        {
            m_lowerList[numSubintervalsIndex] = lower1;
            m_lowerList[maxErrorIndex] = lower2;
            m_upperList[numSubintervalsIndex] = upper1;
            m_integralList[maxErrorIndex] = area2;
            m_integralList[numSubintervalsIndex] = area1;
            list[maxErrorIndex] = error2;
            list[numSubintervalsIndex] = error1;
        }
        else
        {
            m_lowerList[numSubintervalsIndex] = lower2;
            m_upperList[maxErrorIndex] = upper1;
            m_upperList[numSubintervalsIndex] = upper2;
            list[maxErrorIndex] = error1;
            list[numSubintervalsIndex] = error2;
        }

    // Call subroutine quadratureSort to maintain the descending ordering in the list of error estimates and select the
    // subinterval with maxNumberOfIntegrals-th largest error estimate (to be bisected next).
  110  quadratureSort(maxErrorIndex, errorMax, maxNumberOfIntegrals)
    // jump out of do-loop
        if(errorSum <= errorBound) go to 190
    // jump out of do-loop
        if (m_errorCode != 0 || errorSum <= errorBound)
        {
            break;
        }

        if(extrapolationAllowed)
        {
            errorLargest = errorLargest-errorLast
            if(currentLevel+1 <= levmax) errorLargest = errorLargest+erro12
            if(extrapolationPerformed == true) go to 120

    // Test whether the interval to be bisected next is the smallestIntervalLengthest interval.
    if(level(maxErrorIndex)+1 <= levmax) go to 160
    extrapolationPerformed = .true.
    maxNumberOfIntegrals = 2
  120   if(m_errorCode == 3 || errorLargest <= errorTest) go to 140

    // The smallestIntervalLengthest interval has the largest error. Before bisecting decrease the sum of the errors
    // over the larger intervals (errorLargest) and perform extrapolationPerformed.
        id = maxNumberOfIntegrals
        upperBound = m_numSubintervals
        if(m_numSubintervals > (2+m_maxSubintervals/2)) upperBound = m_maxSubintervals+3-m_numSubintervals
        do 130 k = id,upperBound
          maxErrorIndex = listIndices(maxNumberOfIntegrals)
          errorMax = list(maxErrorIndex)
    // jump out of do-loop
          if(level(maxErrorIndex)+1 <= levmax) go to 160
          maxNumberOfIntegrals = maxNumberOfIntegrals+1
  130   continue

    // Perform extrapolationPerformed.
  140   numberElementsList2 = numberElementsList2+1
        m_integralList2(numberElementsList2) = area
        if(numberElementsList2 <= 2) go to 155
        
        qelg(numberElementsList2,m_integralList2,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)
        
        ktmin = ktmin+1
        if(ktmin > 5 && m_estimatedError < 0.001*errorSum) m_errorCode = 5
        if(absoluteEpsilon >= m_estimatedError) go to 150
        ktmin = 0
        m_estimatedError = absoluteEpsilon
        integral = resultEpsilon
        errorCorrection = errorLargest
        errorTest = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(resultEpsilon))
    // jump out of do-loop
        if(m_estimatedError < errorTest) go to 170

    // Prepare bisection of the smallestIntervalLengthest interval.
  150   if(numberElementsList2 == 1) extrapolationAllowed = .true.
        if(m_errorCode >= 5) go to 170
  155   maxErrorIndex = listIndices(1)
        errorMax = list(maxErrorIndex)
        maxNumberOfIntegrals = 1
        extrapolationPerformed = .false.
        levmax = levmax+1
        errorLargest = errorSum
  160 continue

    }// end if extrapolationAllowed
    
    // Set the final integral.
  170 if(m_estimatedError == (std::numeric_limits<Scalar>::max)()) go to 190
      if((m_errorCode+m_errorCode) == 0) go to 180
      if(m_errorCode == 3) m_estimatedError = m_estimatedError+errorCorrection
      if(m_errorCode == 0) m_errorCode = 3
      if(integral != 0. && area != 0.)go to 175
      if(m_estimatedError > errorSum)go to 190
      if(area == 0.) go to 210
      go to 180
  175 if(m_estimatedError/abs(integral) > errorSum/abs(area))go to 190

    // Test on divergence.
  180 if(keySign == (-1) && (std::max)(abs(integral),abs(area)) <= 
     *  absIntegral*0.01) go to 210
      if(0.01 > (integral/area) || (integral/area) > 100. || 
     *  errorSum > abs(area)) m_errorCode = 6
      go to 210

    // Compute global integeral sum.
  190 integral = 0.
      do 200 k = 1,m_numSubintervals
        integral = integral+m_integralList(k)
  200 continue
      m_estimatedError = errorSum
  210 if(m_errorCode > 2) m_errorCode = m_errorCode - 1
      integral = integral*sign
 999  return
      end
