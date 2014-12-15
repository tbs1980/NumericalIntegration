/**
 * \file
 * \brief - This routine calculates an approximation to a given integeral over infinite-bounded intervals
 *          i = integeral of f over (finiteBound,+infiniteBoundsKey), or
 *          i = integeral of f over (-infiniteBoundsKey,finiteBound), or
 *          i = integeral of f over (-infiniteBoundsKey,+infiniteBoundsKey),
 *          hopefully satisfying following claim for accuracy abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i))
 *
 * \details This is a general-purpose, extrapolative, globally adaptive quadrature function for automatic integeration of integrands over infinte intervals.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param finiteBound - Finite bound of integeration range (has no meaning if interval is doubly-infinite-bounded)
 * \param infiniteBoundsKey - Indicating the kind of integeration range involved 
 *                            infiniteBoundsKey = 1 corresponds to (finiteBound,+infiniteBoundsKey),
 *                            infiniteBoundsKey = -1 to  (-infiniteBoundsKey,finiteBound), 
 *                            infiniteBoundsKey = 2 to (-infiniteBoundsKey,+infiniteBoundsKey).
 * \param errorLast - error on the interval currently subdivided (before that subdivision has taken place)
 * \param area - sum of the integerals over the subintervals
 * \param numberOfExtrapolations - number of calls to the extrapolationPerformed routine
 * \param numberElementsList2 - Number of elements currently in m_integralList2. if an appropriate approximation to the compounded integeral has been obtained, it is put in m_integralList2(numberElementsList2) after numberElementsList2 has been increased by one.
 * \param smallestIntervalLength - Length of the smallestIntervalLengthest interval considered up to now, multiplied by 1.5
 * \param errorLargest - sum of the errors over the intervals larger than the smallestIntervalLengthest interval considered up to now
 * \param extrapolationPerformed - logical variable denoting that the routine is attempting to perform extrapolationPerformed. i.e. before subdividing the smallestIntervalLengthest interval we try to decrease the value of errorLargest.
 * \param extrapolationAllowed - logical variable denoting that extrapolationPerformed is no longer allowed (true-value)
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForInfiniteBoundedIntervals(
        const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
        const Scalar desireabsoluteError = Scalar(0.), const Scalar desiredRelativeError = Scalar(0.),
        const QuadratureRule quadratureRule = 1)
{
    // Call error handler if necessary.
    if ((desireabsoluteError <= Scalar(0.) && desiredRelativeError < Eigen::NumTraits<Scalar>::epsilon())
        || m_maxSubintervals < 1 || lenw < m_maxSubintervals * 4)
    {
        m_errorCode = 6;
        return Scalar(0.);
    }

    // Variable initialization
    m_errorCode = 0;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    m_lowerList[0] = lowerLimit;
    m_upperList[0] = upperLimit;
    m_integeralList[0] = Scalar(0.);
    m_errorList[0] = Scalar(0.);
    m_errorListIndices[0] = 0;
    m_errorListIndices[1] = 1;

    Scalar integral = 0.;
    Scalar absDiff = 0.;
    Scalar m_estimatedError = 0.;
    Scalar absIntegral = 0.;


/**************************************** Begin Code Is Different Here *******************************************/

    // First approximation to the integeral.  Determine the interval to be mapped onto (0,1). 
    // If infiniteBoundsKey = 2 the integeral is computed as i = i1+i2, where
    // i1 = integeral of f over (-infiniteBoundsKey,0) and i2 = integeral of f over (0,+infiniteBoundsKey).

    if(infiniteBoundsKey == 2)
    {
        finiteBound = 0.;
    }

    qk15i(f,finiteBound,infiniteBoundsKey,0.,1.,integral,m_estimatedError,absDiff,absIntegral);

    m_numSubintervals = 1;
    m_integralList[0] = integral;
    m_integralList2[1] = integral;
    m_errorList[0] = m_estimatedError;
    m_errorListIndices[1] = 1;

    ktmin = 0;
    numberElementsList2 = 2;

    m_numberOfExtrapolations = 0;
    bool extrapolationPerformed = false;
    bool extrapolationAllowed = false;

    dResult = abs(integral);
    keySign = -1;

    if(dResult >= 1. - (5.*Eigen::NumTraits<Scalar>::epsilon()) * absDiff)
    {
        keySign = 1;
    }

/**************************************** End Code Is Different Here *******************************************/

    // Test on accuracy
    using std::abs;
    Scalar errorBound = (std::max)(desiredAbsoluteError,desiredRelativeError*dResult)

    if (m_maxSubintervals == 1)
    {
        m_errorCode = 1;
    }

    else if (m_estimatedError <= Eigen::NumTraits<Scalar>::epsilon() * Scalar(100.) * absDiff
        && m_estimatedError > errorBound) m_errorCode = 2
    {
        m_errorCode = 2;
    }

    if (m_errorCode != 0
        || (m_estimatedError <= errorBound && m_estimatedError != absIntegral)
        || m_estimatedError == Scalar(0.)) go to 130
    {
        if (quadratureRule == GaussKronrod15)
        {
            m_numEvaluations = m_numEvaluations * 30 - 15;
        }
        else
        {
            m_numEvaluations = (quadratureRule * 10 + 1) * (m_numEvaluations * 2 + 1);
        }
/**************************************** Begin Code Is Different Here *******************************************/
        if(infiniteBoundsKey == 2)
        {
            m_numEvaluations *= 2;
        }
/**************************************** End Code Is Different Here *******************************************/
        return integral;
    }

    // The sum of the integrals over the subintervals.
    Scalar area = integral;

    Scalar errorSum = m_estimatedError;

    // The maximum interval error.
    Scalar errorMax = m_estimatedError;

    // An index into m_errorList at the interval with largest error estimate.
    int maxErrorIndex = 0;
    int maxNumberOfIntegrals = 0;
    int roundOff1 = 0;
    int roundOff2 = 0;
    int roundOff3 = 0;  // This is an additional index over quadratureAdaptive code

    // Main loop for the integration
    for (m_numSubintervals = 2; m_numSubintervals <= m_maxSubintervals; ++m_numSubintervals)
    {
        const int numSubintervalsIndex = m_numSubintervals - 1;
    
        // Bisect the subinterval with the largest error estimate.
        Scalar lower1 = m_lowerList[maxErrorIndex];
        Scalar upper2 = m_upperList[maxErrorIndex];

        Scalar upper1 = (lower1 + upper2) * Scalar(.5);
        Scalar lower2 = upper1;

/**************************************** Begin Code Is Different Here *******************************************/
        errorLast = errorMax;
        call qk15i(f,finiteBound,infiniteBoundsKey,lower1,upper1,area1,error1,absIntegral,absDiff1)
        call qk15i(f,finiteBound,infiniteBoundsKey,lower2,upper2,area2,error2,absIntegral,absDiff2)
/**************************************** End Code Is Different Here *******************************************/

        // Improve previous approximations to integral and error and test for accuracy.
        ++(m_numEvaluations);
        const Scalar area12 = area1 + area2;
        const Scalar error12 = error1 + error2;
        errorSum += error12 - errorMax;
        area += area12 - m_integralList[maxErrorIndex];

        if (absDiff1 != error1 && absDiff2 != error2)
        {
/**************************************** Begin Code Is Different Here *******************************************/
            if(abs(m_integralList[maxErrorIndex] - area12) <= abs(area12) * Scalar(1.e-5)
              && error12 >= errorMax * Scalar(.99) || extrapolationPerformed == false) 
            {
                ++roundOff1;
            }
            if (m_numSubintervals > 10 && error12 > errorMax || extrapolationPerformed == true)
            {
                ++roundOff2;
            }
            if (abs(m_integralList[maxErrorIndex] - area12) <= abs(area12) * Scalar(1.e-5)
              && error12 >= errorMax * Scalar(.99))
            {
                ++roundOff3;
            }
        }
/**************************************** End Code Is Different Here *******************************************/

        m_integralList[maxErrorIndex] = area1;
        m_integralList[numSubintervalsIndex] = area2;

        errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * abs(area));

        if (errorSum > errorBound)
        {
            // Test if the number of subintervals has reached the max allowable.
            if (m_numSubintervals == m_maxSubintervals)
            {
                m_errorCode = 1;
            }
            // Test for roundoff error.
            else if (roundOff1 >= 6 || roundOff2 >= 20)
            {
                m_errorCode = 2;
            }
            // Test for poor integrand behaviour within the integration range.
            else if ((std::max)(abs(lower1), abs(upper2))
                <= (Eigen::NumTraits<Scalar>::epsilon() * Scalar(100.) + Scalar(1.))
                * (abs(lower2) + (std::numeric_limits<Scalar>::min)() * Scalar(1.e3) ))
            {
                m_errorCode = 3;
            }
            // Test for roundoff error during extrapolation
            else if(roundOff2 >= 5)
            {
                m_errorCode = 4;
            }
        } 
          
        // Append the newly-created intervals to the list.
        if (error2 > error1)
        {
            m_lowerList[numSubintervalsIndex] = lower1;
            m_lowerList[maxErrorIndex] = lower2;
            m_upperList[numSubintervalsIndex] = upper1;
            m_integralList[maxErrorIndex] = area2;
            m_integralList[numSubintervalsIndex] = area1;
            m_errorList[maxErrorIndex] = error2;
            m_errorList[numSubintervalsIndex] = error1;
        }
        else
        {
            m_lowerList[numSubintervalsIndex] = lower2;
            m_upperList[maxErrorIndex] = upper1;
            m_upperList[numSubintervalsIndex] = upper2;
            m_errorList[maxErrorIndex] = error1;
            m_errorList[numSubintervalsIndex] = error2;
        }

        // Maintain the descending ordering in the list of error estimates and select the subinterval
        // with the largest error estimate, (the next subinterval to be bisected).
        quadratureSort(maxErrorIndex, errorMax,maxNumberOfIntegrals);
        
        if (m_errorCode != 0 || errorSum <= errorBound)
        {
            break;
        }

/**************************************** Begin Code Is Different Here *******************************************/
        if(m_numSubintervals == 2)
        {
            smallestIntervalLength = 0.375;
            errorLargest = errorSum;
            errorTest = errorBound;
            m_integralList2[2] = area;
        }

        if(extrapolationAllowed == false)
        {
            errorLargest = errorLargest-errorLast
            if(abs(b1-a1) > smallestIntervalLength)
            {
                errorLargest = errorLargest+erro12;
            }
            if(extrapolationPerformed == false)
            {
                // Test whether the interval to be bisected next is the smallestIntervalLengthest interval.
                if(abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex)) > smallestIntervalLength)
                {
                    extrapolationPerformed = true;
                    maxNumberOfIntegrals = 2;
                }

                if(m_errorCode != 3 && errorLargest > errorTest)
                {
                    // The smallestIntervalLengthest interval has the largest error.
                    // Before bisecting, decrease the sum of the errors over the larger
                    // intervals (errorLargest) and perform extrapolation.
                    id = maxNumberOfIntegrals;
                    upperBound = m_numSubintervals;

                    if(m_numSubintervals > (2+m_maxSubintervals/2))
                    {
                        upperBound = m_maxSubintervals+3-m_numSubintervals;
                    }
                    
                    for (size_t k=id; k<upperBound; ++k)
                    {
                        maxErrorIndex = m_errorListIndices(maxNumberOfIntegrals);
                        errorMax = m_errorList(maxErrorIndex);
                        
                        if(abs(m_upperList(maxErrorIndex)-m_lowerList(maxErrorIndex)) > smallestIntervalLength)
                        {
                            break;
                        }

                        ++maxNumberOfIntegrals;
                    }
                }

                // Perform extrapolation.
                ++numberElementsList2;

                m_integralList2(numberElementsList2) = area;

                // Call linear extrapolation Gaussian quadrature routine 
                qelg(numberElementsList2,m_integralList2,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)
                
                ++ktmin;

                if(ktmin > 5 && m_estimatedError < 0.001*errorSum)
                {
                    m_errorCode = 5;
                }
                
                if(absoluteEpsilon < m_estimatedError)
                {
                    ktmin = 0;
                    m_estimatedError = absoluteEpsilon;
                    integral = resultEpsilon;
                    errorCorrection = errorLargest;
                    errorTest = (std::max)(desiredAbsoluteError,desiredRelativeError * abs(resultEpsilon));
                }

                if(m_estimatedError <= errorTest)
                {
                    break;
                }

                // Prepare bisection of the smallestIntervalLengthest interval.
                if(numberElementsList2 == 1)
                {
                    extrapolationAllowed = true;
                }

                if(m_errorCode == 5)
                {
                    break;
                }

                maxErrorIndex = m_errorListIndices(1);
                errorMax = m_errorList(maxErrorIndex);
                maxNumberOfIntegrals = 1;
                extrapolationPerformed = false;
                smallestIntervalLength = smallestIntervalLength * 0.5;
                errorLargest = errorSum;
            }
        }
    }

    // Set final integral and error estimate.
    if(m_estimatedError != (std::numeric_limits<Scalar>::max)() || (m_errorCode != 0)
        
        if(m_errorCode == 3)
        {
            m_estimatedError = m_estimatedError + errorCorrection;
        }

        if(area == 0.)
        {
            m_numEvaluations = 30 * m_numSubintervals - 15;

            if(infiniteBoundsKey == 2)
            {
                m_numEvaluations *= 2;
            }

            return integral;
        }
    }

    if(m_estimatedError/abs(integral) > errorSum/abs(area) || m_estimatedError > errorSum)
    {
        if(keySign != (-1) || (std::max)(abs(integral),abs(area)) > absDiff * 0.01)
        {
            if(0.01 > (integral/area) || (integral/area) > 100 || errorSum > abs(area))
            {
                m_errorCode = 6;
                m_numEvaluations = 30 * m_numSubintervals - 15;

                if(infiniteBoundsKey == 2)
                {
                    m_numEvaluations *= 2;
                }
                return integral;
            }
        }
    }

    if(infiniteBoundsKey == 2)
    {
        m_numEvaluations *= 2;
    }

/**************************************** End Code Is Different Here *******************************************/
    integral = Scalar(0.);

    for (int k = 0; k < m_numSubintervals; ++k)
    {
        integral += m_integralList[k];
    }

    m_estimatedError = errorSum;

    if (quadratureRule == GaussKronrod15)
    {
        m_numEvaluations = m_numEvaluations * 30 - 15;
    }
    else
    {
        m_numEvaluations = (quadratureRule * 10 + 1) * (m_numEvaluations * 2 + 1);
    }

    return integral;
}
