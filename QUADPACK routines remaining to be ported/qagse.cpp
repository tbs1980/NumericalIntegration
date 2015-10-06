/**
 * \file
 * \brief - This routine calculates an approximation integral to a given definite integeral i = integeral of f over (lowerLimit, upperLimit),
 *          hopefully satisfying following claim for accuracy abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i)).
 *
 * \details - This function is an adaptive quadrature routine for the computation of definite integerals inclusive of simgularities.
 *            It is a general-purpose, globally-adaptive, extrapolative, automatic integerator for functions with end point singularities.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] desiredAbsoluteError - Absolute accuracy requested
 * \param[] desiredRelativeError - Relative accuracy requested. If  desiredAbsoluteError <= 0, and desiredRelativeError < max(50*rel.mach.acc.,0.5d-28), the routine will end with m_errorCode = 6.
 * \param[] integral - Approximation to the integeral
 * \param[] m_estimatedError - Estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)
 * \param[] maxErrorIndex - Pointer to the interval with largest error estimate
 * \param[] errorMax - list(maxErrorIndex)
 * \param[] errorLast - Error on the interval currently subdivided (before that subdivision has taken place)
 * \param[] area - sum of the integerals over the subintervals
 * \param[] errorSum - sum of the errors over the subintervals
 * \param[] errorBound - requested accuracy max(desiredAbsoluteError,desiredRelativeError * abs(integral))
 * \param[] numberOfExtrapolations - number of calls to the extrapolationPerformed routine
 * \param[] numberElementsList2 - number of elements currently in m_integralList2. if an appropriate approximation to the compounded integeral has been obtained it is put in m_integralList2(numberElementsList2) after numberElementsList2 has been increased by one.
 * \param[] smallestIntervalLength - length of the smallestIntervalLengthest interval considered up to now, multiplied by 1.5
 * \param[] errorLargest - sum of the errors over the intervals larger than the smallestIntervalLengthest interval considered up to now
 * \param[] extrapolationPerformed - logical variable denoting that the routine is attempting to perform extrapolationPerformed. i.e. before subdividing the smallestIntervalLengthest interval we try to decrease the value of errorLargest.
 * \param[] extrapolationAllowed - logical variable denoting that extrapolationPerformed is still allowed (true value)
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForEndPointSingularities(
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

    Scalar a;
    Scalar absoluteEpsilon;
    Scalar absIntegral;
    Scalar area;
    Scalar area1;
    Scalar area12;
    Scalar area2;
    Scalar a1;
    Scalar a2;
    Scalar b;
    Scalar b1;
    Scalar b2;
    Scalar errorCorrection;
    Scalar absDiff;
    Scalar absDiff1;
    Scalar absDiff2;
    Scalar dResult;
    Scalar desiredAbsoluteError;
    Scalar desiredRelativeError;
    Scalar errorLargest;
    Scalar errorLast;
    Scalar errorBound;
    Scalar errorMax;
    Scalar error1;
    Scalar error2;
    Scalar erro12;
    Scalar errorSum;
    Scalar errorTest;
    Scalar resultEpsilon;
    Scalar res3la;
    Scalar r1mach;
    Scalar smallestIntervalLength;

    int id;
    int roundOff1;
    int roundOff2;
    int roundOff3;
    int upperBound;
    int k;
    int keySign;
    int ktmin;
    int numberOfExtrapolations;
    int maxNumberOfIntegrals;
    int numberElementsList2;

    res3la(3);

    m_errorCode = 0;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    integral = 0.;
    m_estimatedError = 0.;
      
    if ((desireabsoluteError <= Scalar(0.) && desiredRelativeError < (std::max)(50. * Eigen::NumTraits<Scalar>::epsilon(),0.5e-14) || m_maxSubintervals < 1)
    {
        m_errorCode = 6;
        return Scalar(0.);
    }

    // First approximation to the integeral
    m_errorCode = 0
    qk21(f,lowerLimit, upperLimit,integral,m_estimatedError,absDiff,absIntegral)

    // Test on accuracy.
    dResult = abs(integral)
    errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*dResult)
    m_numSubintervals = 1
    m_integralList(1) = integral
    list(1) = m_estimatedError
    listIndices(1) = 1
    if(m_estimatedError <= 100.*Eigen::NumTraits<Scalar>::epsilon()*absDiff && m_estimatedError > errorBound) m_errorCode = 2
    if(m_maxSubintervals == 1) m_errorCode = 1
    if(m_errorCode != 0 || (m_estimatedError <= errorBound && m_estimatedError != absIntegral) || m_estimatedError == 0.) go to 140

    m_integralList2(1) = integral
    errorMax = m_estimatedError
    maxErrorIndex = 1
    area = integral
    errorSum = m_estimatedError
    m_estimatedError = (std::numeric_limits<Scalar>::max)()
    maxNumberOfIntegrals = 1
    numberOfExtrapolations = 0
    numberElementsList2 = 2
    ktmin = 0;
    bool extrapolationPerformed = false;
    bool extrapolationAllowed = true;
    roundOff1 = 0;
    roundOff2 = 0;
    roundOff3 = 0;

    keySign = -1
    if(dResult >= (1.-5.*Eigen::NumTraits<Scalar>::epsilon())*absDiff) keySign = 1

    // Main loop for the integration
    do 90 m_numSubintervals = 2,m_maxSubintervals

    // Bisect the subinterval with the maxNumberOfIntegrals-th largest error estimate.
    a1 = m_lowerList(maxErrorIndex)
    b1 = 0.5*(m_lowerList(maxErrorIndex)+m_upperList(maxErrorIndex))
    a2 = b1
    b2 = m_upperList(maxErrorIndex)
    errorLast = errorMax
    call qk21(f,a1,b1,area1,error1,absIntegral,absDiff1)
    call qk21(f,a2,b2,area2,error2,absIntegral,absDiff2)

    // Improve previous approximations to integeral and error and test for accuracy.
    area12 = area1+area2
    erro12 = error1+error2
    errorSum = errorSum+erro12-errorMax
    area = area+area12-m_integralList(maxErrorIndex)
    if(absDiff1 == error1 || absDiff2 == error2) go to 15
    if(abs(m_integralList(maxErrorIndex)-area12) > 0.1e-04*abs(area12) || erro12 < 0.99 *errorMax) go to 10
    if(extrapolationPerformed) roundOff2 = roundOff2+1
    if(.not.extrapolationPerformed) roundOff1 = roundOff1+1
10  if(m_numSubintervals > 10 && erro12 > errorMax) roundOff3 = roundOff3+1
15  m_integralList(maxErrorIndex) = area1
    m_integralList(m_numSubintervals) = area2
    errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(area))

    // Test for roundoff error and eventually set error flag.
    if(roundOff1+roundOff2 >= 10 || roundOff3 >= 20) m_errorCode = 2
    if(roundOff2 >= 5) m_errorCode = 3

    // Set error flag in the case that the number of subintervals equals m_maxSubintervals.
    if(m_numSubintervals == m_maxSubintervals) m_errorCode = 1

    // Set error flag in the case of bad integerand behaviour at a point of the integeration range.
    if((std::max)(abs(a1),abs(b2)) <= (1.+100.*Eigen::NumTraits<Scalar>::epsilon())*(abs(a2)+0.1e+04*(std::numeric_limits<Scalar>::min)())) m_errorCode = 4

    // Append the newly-created intervals to the list.
    if(error2 > error1) go to 20
    m_lowerList(m_numSubintervals) = a2
    m_upperList(maxErrorIndex) = b1
    m_upperList(m_numSubintervals) = b2
    list(maxErrorIndex) = error1
    list(m_numSubintervals) = error2
    go to 30
20  m_lowerList(maxErrorIndex) = a2
    m_lowerList(m_numSubintervals) = a1
    m_upperList(m_numSubintervals) = b1
    m_integralList(maxErrorIndex) = area2
    m_integralList(m_numSubintervals) = area1
    list(maxErrorIndex) = error2
    list(m_numSubintervals) = error1

    // Call subroutine quadratureSort to maintain the descending ordering in the list of error estimates and select the 
    // subinterval with maxNumberOfIntegrals-th largest error estimate (to be bisected next).
30  call quadratureSort(m_maxSubintervals,m_numSubintervals,maxErrorIndex,errorMax,list,listIndices,maxNumberOfIntegrals)

    // jump out of do-loop
    if(errorSum <= errorBound) go to 115
    // jump out of do-loop
    if(m_errorCode != 0) go to 100
    if(m_numSubintervals == 2) go to 80
    if(extrapolationAllowed)
    {
    errorLargest = errorLargest-errorLast
    if(abs(b1-a1) > smallestIntervalLength) errorLargest = errorLargest+erro12
    if(extrapolationPerformed) go to 40

    // Test whether the interval to be bisected next is the smallestIntervalLengthest interval.
    if(abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex)) > smallestIntervalLength) go to 90
    extrapolationPerformed = .true.
    maxNumberOfIntegrals = 2
40  if(m_errorCode == 3 || errorLargest <= errorTest) go to 60

    // The smallestIntervalLengthest interval has the largest error. Before bisecting decrease the sum of the errors
    // over the larger intervals (errorLargest) and perform extrapolationPerformed.
    id = maxNumberOfIntegrals
    upperBound = m_numSubintervals
    if(m_numSubintervals > (2+m_maxSubintervals/2)) upperBound = m_maxSubintervals+3-m_numSubintervals
    do 50 k = id,upperBound
      maxErrorIndex = listIndices(maxNumberOfIntegrals)
      errorMax = list(maxErrorIndex)
    // jump out of do-loop
    if(abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex)) > smallestIntervalLength) go to 90
    maxNumberOfIntegrals = maxNumberOfIntegrals+1
50  continue

    // Perform extrapolationPerformed.
60  numberElementsList2 = numberElementsList2+1
    m_integralList2(numberElementsList2) = area
    call qelg(numberElementsList2,m_integralList2,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)
    ktmin = ktmin+1
    if(ktmin > 5 && m_estimatedError < 0.001*errorSum) m_errorCode = 5
    if(absoluteEpsilon >= m_estimatedError) go to 70
    ktmin = 0
    m_estimatedError = absoluteEpsilon
    integral = resultEpsilon
    errorCorrection = errorLargest
    errorTest = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(resultEpsilon))
    
    // jump out of do-loop
    if(m_estimatedError <= errorTest) go to 100

    // Prepare bisection of the smallestIntervalLengthest interval.
70  if(numberElementsList2 == 1) extrapolationAllowed = .true.
    if(m_errorCode == 5) go to 100
    maxErrorIndex = listIndices(1)
    errorMax = list(maxErrorIndex)
    maxNumberOfIntegrals = 1
    extrapolationPerformed = .false.
    smallestIntervalLength = smallestIntervalLength*0.5
    errorLargest = errorSum
    go to 90
80  smallestIntervalLength = abs(upperLimit-lowerLimit)*0.375 
    errorLargest = errorSum
    errorTest = errorBound
    m_integralList2(2) = area
90  continue
    }//end if extrapolation allowed

    // Set final integral and error estimate.
100 if(m_estimatedError == (std::numeric_limits<Scalar>::max)()) go to 115
    if(m_errorCode+m_errorCode == 0) go to 110
    if(m_errorCode == 3) m_estimatedError = m_estimatedError+errorCorrection
    if(m_errorCode == 0) m_errorCode = 3
    if(integral != 0. && area != 0.) go to 105
    if(m_estimatedError > errorSum) go to 115
    if(area == 0.) go to 130
    go to 110
105 if(m_estimatedError/abs(integral) > errorSum/abs(area)) go to 115

    // Test on divergence.
110 if(keySign == (-1) && (std::max)(abs(integral),abs(area)) <= absDiff*0.01) go to 130
      if(0.01 > (integral/area) || (integral/area) > 100. || errorSum > abs(area)) m_errorCode = 6
      go to 130

//   Compute global integeral sum.
115 integral = 0.
    do 120 k = 1,m_numSubintervals
    integral = integral+m_integralList(k)
120 continue
    m_estimatedError = errorSum
130 if(m_errorCode > 2) m_errorCode = m_errorCode-1
140 m_numEvaluations = 42*m_numSubintervals-21
999 return
    end
