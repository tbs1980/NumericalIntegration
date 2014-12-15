/**
 * \file
 * \brief This function calculates an approximation I' to a given definite integeral I, the integeral of f*w from the lowerLimit to upperLimit, (where w shows a singular behaviour
 *        at the end points see parameter quadratureRule), hopefully satisfying abs(I - I') <= max(desireabsoluteError, desiredRelativeError * abs(I)).
 *
 * \details This file contains routines for numerical integeration for inegrands with singularities or discontinuities using Curtis-Clenshaw quadrature.
 *          This routine is better suited than quadratureAdaptive() for integerands with singularities or discontinuities.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[in/out] f The function defining the integerand function.
 * \param[in] lowerLimit The lower m_maxSubintervals of integeration.
 * \param[in] upperLimit The upper m_maxSubintervals of integeration.
 * \param[in] desireabsoluteError The absolute accuracy requested.
 * \param[in] desiredRelativeError The relative accuracy requested.
 *            If desireabsoluteError <= 0 and desiredRelativeError < 50 * machinePrecision,
 *            the routine will end with errorCode = 6.
 * \param[in] quadratureRule The local Gauss-Kronrod quadrature rule to use.
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForSingularities(
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
    m_lowerList[0] = lowerLimit;
    m_upperList[0] = upperLimit;
    m_integralList[0] = Scalar(0.);
    m_errorList[0] = Scalar(0.);
    m_errorListIndices[0] = 0;
    m_errorListIndices[1] = 1;

    Scalar integral = 0.;
    Scalar absDiff = 0.;
    Scalar m_estimatedError = 0.;
    Scalar absIntegral = 0.;


/**************************************** Begin Code Is Different Here *******************************************/
    Scalar alpha = 0.;
    Scalar beta = 0.;

    // Compute the modified chebyshev moments.
    dqmomo(alphlowerLimit, beta,ri,rj,rg,rh,quadratureRule)

    // Half-length and midpoint of the interval.
    const Scalar halfLength = (upperLimit - lowerLimit) * Scalar(.5);
    const Scalar center = (lowerLimit + upperLimit) * Scalar(.5);
    
    const Scalar  fCenter = f(center);

    // First approximation to the integeral by integerating over the half-intervals.
    dqc25s(f,lowerLimit,upperLimit,lowerLimit,center,alphlowerLimit, beta,ri,rj,rg,rh,area1,error1,approximateIntegralResults1,quadratureRule,nEval)
    ++(m_numEvaluations);

    dqc25s(f,lowerLimit,upperLimit,center,upperLimit,alphlowerLimit, beta,ri,rj,rg,rh,area2,error2,approximateIntegralResults2,quadratureRule,nEval)
    ++(m_numEvaluations);

    m_numSubintervals = 2;
    area12 = area1 + area2;
    m_estimatedError = error1+error2;

    // integeral and error list initializations
    if(error1 <= error2)
    {
        lowerintegeralList(0) = lowerLimit
        lowerintegeralList(1) = center
        upperintegeralList(0) = center
        upperintegeralList(1) = upperLimit
        integeralList(0) = area1
        integeralList(1) = area2
        errorList(0) = error1
        errorList(1) = error2
    }
    else
    {
        lowerintegeralList(0) = center
        lowerintegeralList(1) = lowerLimit
        upperintegeralList(0) = upperLimit
        upperintegeralList(1) = center
        integeralList(0) = area2
        integeralList(1) = area1
        errorList(0) = error2
        errorList(1) = error1
    }

/**************************************** End Code Is Different Here *******************************************/

    m_numSubintervals = 1;
    m_integralList[0] = integral;
    m_errorList[0] = m_estimatedError;

    // Test on accuracy.
    using std::abs;
    Scalar errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * abs(integral));

    if (m_maxSubintervals == 1)
    {
        m_errorCode = 1;
    }
    else if (m_estimatedError <= Eigen::NumTraits<Scalar>::epsilon() * Scalar(50.) * absDiff
        && m_estimatedError > errorBound)
    {
        m_errorCode = 2;
    }

    if (m_errorCode != 0
        || (m_estimatedError <= errorBound && m_estimatedError != absIntegral)
        || m_estimatedError == Scalar(0.))
    {
        if (quadratureRule == GaussKronrod15)
        {
            m_numEvaluations = m_numEvaluations * 30 + 15;
        }
        else
        {
            m_numEvaluations = (quadratureRule * 10 + 1) * (m_numEvaluations * 2 + 1);
        }

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

    // Main loop for the integration
    for (m_numSubintervals = 2; m_numSubintervals <= m_maxSubintervals; ++m_numSubintervals)
    {
        const int numSubintervalsIndex = m_numSubintervals - 1;

        // Bisect the subinterval with the largest error estimate.
        Scalar lower1 = m_lowerList[maxErrorIndex];
        Scalar upper2 = m_upperList[maxErrorIndex];

        Scalar upper1 = (lower1 + upper2) * Scalar(.5);
        Scalar lower2 = upper1;

        Scalar error1;
        Scalar error2;
        Scalar absDiff1;
        Scalar absDiff2;


/**************************************** Begin Code Is Different Here *******************************************/
        call dqc25s(f,lowerLimit,upperLimit,lowerLimit,center,alphlowerLimit, beta,ri,rj,rg,rh,area1,error1,approximateIntegralResults1,quadratureRule,nEval)
        ++(m_numEvaluations);

        call dqc25s(f,lowerLimit,upperLimit,center,upperLimit,alphlowerLimit, beta,ri,rj,rg,rh,area2,error2,approximateIntegralResults2,quadratureRule,nEval)
        ++(m_numEvaluations);
/**************************************** End Code Is Different Here *******************************************/

        // Improve previous approximations to integral and error and test for accuracy.
        ++(m_numEvaluations);
        const Scalar area12 = area1 + area2;
        const Scalar error12 = error1 + error2;
        errorSum += error12 - errorMax;
        area += area12 - m_integralList[maxErrorIndex];

        if (absDiff1 != error1 && absDiff2 != error2)
        {
            if (abs(m_integralList[maxErrorIndex] - area12) <= abs(area12) * Scalar(1.e-5)
              && error12 >= errorMax * Scalar(.99))
            {
                ++roundOff1;
            }
            if (m_numSubintervals > 10 && error12 > errorMax)
            {
                ++roundOff2;
            }
        }

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
        quadratureSort(maxErrorIndex, errorMax, maxNumberOfIntegrals);

        if (m_errorCode != 0 || errorSum <= errorBound)
        {
            break;
        }
    }

    integral = Scalar(0.);

    for (int k = 0; k < m_numSubintervals; ++k)
    {
        integral += m_integralList[k];
    }

    m_estimatedError = errorSum;

    if (quadratureRule == GaussKronrod15)
    {
        m_numEvaluations = m_numEvaluations * 30 + 15;
    }
    else
    {
        m_numEvaluations = (quadratureRule * 10 + 1) * (m_numEvaluations * 2 + 1);
    }

    return integral;
}