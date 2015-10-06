/**
 * \file
 * \brief - The routine calculates an approximation integral to a given definite integeral i = integeral of f(x)*w(x) 
 *            over (lowerLimit, upperLimit) where w(x) = cos(omega*x) or w(x) = sin(omega*x), hopefully satisfying the
 *            claim for accuracy abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i)).
 *
 * \details - This routine is a globally adaptive, extrapolative, Clenshaw-curtis method, automatic,  integerator, special-purpose, integerand with cos or sin factor oscillatory integerals with (end point) singularities.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] omega  - Parameter in the integerand weight function
 * \param[] integer - Indicates which of the weight functions is to be used
 *             integer = 1      w(x) = cos(omega*x)
 *             integer = 2      w(x) = sin(omega*x)
 *             if integer != 1 and integer != 2, the routine will end with m_errorCode = 6.
 * \param[] desiredAbsoluteError - Absolute accuracy requested
 * \param[] desiredRelativeError - Relative accuracy requested.  If  desiredAbsoluteError <= 0 and desiredRelativeError < max(50*rel.mach.acc.,0.5d-28), the routine will end with m_errorCode = 6.
 * \param[] icall - If dqawoe is to be used only once, icall must be set to 1.  assume that during this call, the Chebyshev moments (for clenshaw-curtis integeration of degree 24) have been computed for intervals of
 *                     lenghts (abs(upperLimit-lowerLimit))*2**(-l), l=0,1,2,...momentsComputed-1. If icall > 1 this means that dqawoe has been called twice or more on intervals of the same
 *                     length abs(upperLimit-lowerLimit). the chebyshev moments already computed are then re-used in subsequent calls.  If icall < 1, the routine will end with m_errorCode = 6.
 * \param[] maxp1 - Gives an upper finiteBound on the number of Chebyshev moments which can be stored, i.e. for the intervals of lenghts abs(upperLimit-lowerLimit)*2**(-l), from l=0 to maxp1-2, maxp1 >= 1. If maxp1 < 1, the routine will end with m_errorCode = 6.
 * \param[] integral - Approximation to the integeral
 * \param[] m_estimatedError - Estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)
 * \param[] nnlog - Vector of dimension at least m_maxSubintervals, containing the subdivision levels of the subintervals, i.e. iwork(i) = l means that the subinterval numbered i is of length abs(upperLimit-lowerLimit)*2**(1-l)
 * \param[] momentsComputed - Indicating that the chebyshev moments have been computed for intervals of lengths (abs(upperLimit-lowerLimit))*2**(-l), from l=0 to momentsComputed-1, momentsComputed < maxp1
 * \param[] chebyshevMoments - Array of dimension (maxp1,25) containing the Chebyshev moments.  The dimension of m_integralList2 is determined by  the value of limexp in subroutine qelg (m_integralList2 should be of dimension (limexp+2) at least).
 * \param[] maxErrorIndex - pointer to the interval with largest error estimate
 * \param[] errorMax - list(maxErrorIndex)
 * \param[] errorLast - error on the interval currently subdivided
 * \param[] area - sum of the integerals over the subintervals
 * \param[] errorSum - sum of the errors over the subintervals
 * \param[] errorBound - requested accuracy max(desiredAbsoluteError,desiredRelativeError*abs(integral))
 * \param[] numberOfExtrapolations - number of calls to the extrapolationPerformed routine
 * \param[] numberElementsList2 - number of elements in m_integralList2. if an appropriate approximation to the compounded integeral has been obtained it is put in m_integralList2(numberElementsList2) after numberElementsList2 has been increased by one
 * \param[] smallestIntervalLength - length of the smallestIntervalLengthest interval considered up to now, multiplied by 1.5
 * \param[] errorLargest - sum of the errors over the intervals larger than the smallestIntervalLengthest interval considered up to now
 * \param[] extrapolationPerformed - logical variable denoting that the routine is attempting to perform extrapolationPerformed, i.e. before subdividing the smallestIntervalLengthest interval we try to decrease the value of errorLargest
 * \param[] extrapolationAllowed - logical variable denoting that extrapolationPerformed is still allowed (true value)
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForOscillatoryBehaviorWithEndPointSingularities(
        const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
        const Scalar desireabsoluteError = Scalar(0.), const Scalar desiredRelativeError = Scalar(0.),
        const QuadratureRule quadratureRule = 1)
{
    if ((desiredAbsoluteError <= Scalar(0.) && desiredRelativeError < Eigen::NumTraits<Scalar>::epsilon())
        || m_maxSubintervals < 1
        || (integer != 1 && integer != 2)
        || icall < 1
        || maxp1 < 1)
    {
        m_errorCode = 6;
        return Scalar(0.);
    }

    m_errorCode = 0;
    m_numEvaluations = 0;
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

      Scalar absoluteEpsilon;
      Scalar area;
      Scalar area1;
      Scalar integral;
      Scalar area2;
      Scalar a1;
      Scalar a2;
      Scalar b1;
      Scalar b2;
      Scalar chebyshevMoments;
      Scalar errorCorrection;
      Scalar absDiff1;
      Scalar absDiff2;
      Scalar absDiff;
      Scalar omega;
      Scalar dResult;
      Scalar desiredAbsoluteError;
      Scalar desiredRelativeError;
      Scalar errorLargest;
      Scalar errorLast;
      Scalar errorBound;
      Scalar errorMax;
      Scalar errorTest;
      Scalar omega;
      Scalar absIntegral;
      Scalar resultEpsilon;
      Scalar integral;
      Scalar res3la;
      Scalar smallestIntervalLength;
      Scalar width;
      
      int icall;
      int id;
      int integer;
      int roundOff1;
      int roundOff2;
      int roundOff3;
      int upperBound;
      int keySign;
      int ktmin;
      int maxErrorIndex;
      int maxp1;
      int momentsComputed;
      int nEval;

      int numberOfExtrapolations;
      int maxNumberOfIntegrals;
      int nrMoments;
      int numberElementsList2;

      bool extrapolationPerformed;
      bool extrapolationAllowed;
      bool extAll;

      res3la[3]
      chebyshevMoments[maxp1,25]
      nnlog[m_maxSubintervals];

      int nnlog;
      nnlog(1) = 0;
    //first approximation to the integeral
      omega = abs(omega)
      nrMoments = 0
    if (icall <= 1)
    {
        momentsComputed = 0;
    }

    qc25f(f,lowerLimit, upperLimit,omega,integer,nrMoments,maxp1,0,integral,m_estimatedError,m_numEvaluations,absDiff,absIntegral,momentsComputed,chebyshevMoments)

    // Test on accuracy.
      dResult = abs(integral);
      errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*dResult)
      m_integralList(1) = integral
      list(1) = m_estimatedError
      listIndices(1) = 1
      if(m_estimatedError <= 100.*Eigen::NumTraits<Scalar>::epsilon()*absDiff && m_estimatedError > errorBound) m_errorCode = 2
      if(m_maxSubintervals == 1) m_errorCode = 1
      if(m_errorCode != 0 || m_estimatedError <= errorBound) go to 200

    // Initializations
    (std::numeric_limits<Scalar>::min)() = r1mach(1)
    (std::numeric_limits<Scalar>::max)() = r1mach(2)
    errorMax = m_estimatedError
    maxErrorIndex = 1
    area = integral
    Scalar errorSum = m_estimatedError;
    m_estimatedError = (std::numeric_limits<Scalar>::max)()
    maxNumberOfIntegrals = 1
    extrapolationPerformed = false;
    extrapolationAllowed = true;
    m_errorCode = 0
    roundOff1 = 0
    roundOff2 = 0
    roundOff3 = 0
    ktmin = 0
    smallestIntervalLength = abs(upperLimit-lowerLimit)*0.75 
    numberOfExtrapolations = 0
    numberElementsList2 = 0
    extAll = .false.
    if(0.5*abs(upperLimit-lowerLimit)*omega > 2.) go to 10
    numberElementsList2 = 1
    extAll = .true.
    m_integralList2(1) = integral
    10 if(0.25 *abs(upperLimit-lowerLimit)*omega <= 2.) extAll = .true.
    keySign = -1
    if(dResult >= (1.-5.*Eigen::NumTraits<Scalar>::epsilon())*absDiff) keySign = 1

    // mMain loop
      do 140 m_numSubintervals = 2,m_maxSubintervals

    // Bisect the subinterval with the maxNumberOfIntegrals-th largest error estimate.
        nrMoments = nnlog(maxErrorIndex)+1
        a1 = m_lowerList(maxErrorIndex)
        b1 = 0.5*(m_lowerList(maxErrorIndex)+m_upperList(maxErrorIndex))
        a2 = b1
        b2 = m_upperList(maxErrorIndex)
        errorLast = errorMax
        
        Scalar error1;
        Scalar error2;

        qc25f(f,a1,b1,omega,integer,nrMoments,maxp1,0,area1,error1,nEval,absIntegral,absDiff1,momentsComputed,chebyshevMoments)
        m_numEvaluations = m_numEvaluations+nEval
        
        qc25f(f,a2,b2,omega,integer,nrMoments,maxp1,1,area2,error2,nEval,absIntegral,absDiff2,momentsComputed,chebyshevMoments)
        m_numEvaluations = m_numEvaluations+nEval

    // Improve previous approximations to integeral and error and test for accuracy.
        Scalar area12 = area1 + area2;
        Scalar erro12 = error1 + error2;
        errorSum += erro12 - errorMax;
        area += integral - m_integralList(maxErrorIndex);

        if(absDiff1 != error1 && absDiff2 != error2)
        {
            if(abs(m_integralList(maxErrorIndex)-integral) <= 0.00001*abs(integral) && erro12 >= 0.99 *errorMax)
            {
                if(extrapolationPerformed)
                {
                    ++roundOff2;
                }
                else
                {
                    ++roundOff1;
                }
            }
        }

   20   if(m_numSubintervals > 10 && erro12 > errorMax) roundOff3 = roundOff3+1
   25   m_integralList(maxErrorIndex) = area1
        m_integralList(m_numSubintervals) = area2
        nnlog(maxErrorIndex) = nrMoments
        nnlog(m_numSubintervals) = nrMoments
        errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(area))

    // Test for roundoff error and eventually set error flag
        if(roundOff1+roundOff2 >= 10 || roundOff3 >= 20) m_errorCode = 2
        if(roundOff2 >= 5) m_errorCode = 3

    // Set error flag in the case that the number of subintervals equals m_maxSubintervals.
        if(m_numSubintervals == m_maxSubintervals) m_errorCode = 1

    // Set error flag in the case of bad integerand behaviour at a point of the integeration range.
    if((std::max)(abs(a1),abs(b2)) <= (1.+100.*Eigen::NumTraits<Scalar>::epsilon())*(abs(a2)+0.1e+04*(std::numeric_limits<Scalar>::min)())) m_errorCode = 4

    // Append the newly-created intervals to the list.
        if(error2 > error1) go to 30
        m_lowerList(m_numSubintervals) = a2
        m_upperList(maxErrorIndex) = b1
        m_upperList(m_numSubintervals) = b2
        list(maxErrorIndex) = error1
        list(m_numSubintervals) = error2
        go to 40
   30   m_lowerList(maxErrorIndex) = a2
        m_lowerList(m_numSubintervals) = a1
        m_upperList(m_numSubintervals) = b1
        m_integralList(maxErrorIndex) = area2
        m_integralList(m_numSubintervals) = area1
        list(maxErrorIndex) = error2
        list(m_numSubintervals) = error1

    // Call subroutine quadratureSort to maintain the descending ordering in the list of error estimates and select
    // the subinterval with maxNumberOfIntegrals-th largest error estimate (to be bisected next).
   40   call quadratureSort(m_maxSubintervals,m_numSubintervals,maxErrorIndex,errorMax,list,listIndices,maxNumberOfIntegrals)
    // jump out of do-loop
      if(errorSum <= errorBound) go to 170
      if(m_errorCode != 0) go to 150
        if(m_numSubintervals == 2 && extAll) go to 120
        if(extrapolationAllowed)
        {
        if(.not.extAll) go to 50
        errorLargest = errorLargest-errorLast
        if(abs(b1-a1) > smallestIntervalLength) errorLargest = errorLargest+erro12
        if(extrapolationPerformed) go to 70

    // Test whether the interval to be bisected next is the smallestIntervalLengthest interval.
   50   width = abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex))
        if(width > smallestIntervalLength) go to 140
        if(extAll) go to 60

    // Test whether we can start with the extrapolationPerformed procedure (we do this if we integerate over the
    // next interval with use of a gauss-kronrod rule - see subroutine qc25f).
        smallestIntervalLength = smallestIntervalLength*0.5
        if(0.25 *width*omega > 2.) go to 140
        extAll = .true.
        go to 130
   60   extrapolationPerformed = .true.
        maxNumberOfIntegrals = 2
   70   if(m_errorCode == 3 || errorLargest <= errorTest) go to 90

    // The smallestIntervalLengthest interval has the largest error. Before bisecting decrease the sum of the errors
    // over the larger intervals (errorLargest) and perform extrapolationPerformed.
        upperBound = m_numSubintervals
        if (m_numSubintervals > (m_maxSubintervals/2+2)) upperBound = m_maxSubintervals+3-m_numSubintervals
        id = maxNumberOfIntegrals
        do 80 k = id,upperBound
          maxErrorIndex = listIndices(maxNumberOfIntegrals)
          errorMax = list(maxErrorIndex)
          if(abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex)) > smallestIntervalLength) go to 140
          maxNumberOfIntegrals = maxNumberOfIntegrals+1
   80   continue

    // Perform extrapolationPerformed.
   90   numberElementsList2 = numberElementsList2+1
        m_integralList2(numberElementsList2) = area
        if(numberElementsList2 < 3) go to 110
        call qelg(numberElementsList2,m_integralList2,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)
        ktmin = ktmin+1
        if(ktmin > 5 && m_estimatedError < 0.001*errorSum) m_errorCode = 5
        if(absoluteEpsilon >= m_estimatedError) go to 100
        ktmin = 0
        m_estimatedError = absoluteEpsilon
        integral = resultEpsilon
        errorCorrection = errorLargest
        errorTest = (std::max)(desiredAbsoluteError,desiredRelativeError*abs(resultEpsilon))
    // Jump out of do-loop
        if(m_estimatedError <= errorTest) go to 150

    // Prepare bisection of the smallestIntervalLengthest interval.
  100   if(numberElementsList2 == 1) extrapolationAllowed = false;
        if(m_errorCode == 5) go to 150
  110   maxErrorIndex = listIndices(1)
        errorMax = list(maxErrorIndex)
        maxNumberOfIntegrals = 1
        extrapolationPerformed = .false.
        smallestIntervalLength = smallestIntervalLength*0.5
        errorLargest = errorSum
        go to 140
  120   smallestIntervalLength = smallestIntervalLength*0.5
        numberElementsList2 = numberElementsList2+1
        m_integralList2(numberElementsList2) = area
  130   errorTest = errorBound
        errorLargest = errorSum
  140 continue
    }// end if extrapolationAllowed

    // Set the final integral.
  150 if(m_estimatedError == (std::numeric_limits<Scalar>::max)() || numberOfExtrapolations == 0) go to 170
      if(m_errorCode+m_errorCode == 0) go to 165
      if(m_errorCode == 3) m_estimatedError = m_estimatedError+errorCorrection
      if(m_errorCode == 0) m_errorCode = 3
      if(integral != 0. && area != 0.) go to 160
      if(m_estimatedError > errorSum) go to 170
      if(area == 0.) go to 190
      go to 165
  160 if(m_estimatedError/abs(integral) > errorSum/abs(area)) go to 170

    // Test on divergence.
165 if(keySign == (-1) && (std::max)(abs(integral),abs(area)) <= absDiff*0.01)
    {
    
        if(0.01 > (integral/area) || (integral/area) > 100. || errorSum >= abs(area))
        {
            m_errorCode = 5;
        }
        
        if (integer == 2 && omega < 0.) 
        {
            integral *= -1.;
        }

        return integral;
    }

    // Compute global integeral sum.
170 integral = 0.;

    for (Size_t k=0; k<m_numSubintervals; ++k)
    {
        integral = integral+m_integralList(k);
    }
180   
    m_estimatedError = errorSum;

190 if (m_errorCode > 2)
    {
        --m_errorCode;
    }
200 
if (integer == 2 && omega < 0.) 
    {
        integral *= -1.;
    }
  999 return
      end
