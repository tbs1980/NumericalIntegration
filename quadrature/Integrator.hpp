/**
 * \file Integrator.h
 * \sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package
 *     for Automatic Integration, Springer Verlag, 1983.
 */

#ifndef EIGEN_INTEGRATOR_H
#define EIGEN_INTEGRATOR_H

namespace Eigen
{


/**
 * \ingroup Quadrature_Module
 *
 * \brief This class provides numerical integration functionality.
 *
 * Memory management and additional information (e.g. number of function calls, estimated error),
 * are provided by the class in a way that can support the future porting of additional quadrature
 * functions from the QUADPACK library.
 *
 * \todo Ensure only appropriates types are used for Scalar_, e.g. prohibit integers.
 */
template <typename Scalar_>
class Integrator
{
public:
  typedef Scalar_ Scalar;

  /**
   * \brief The local Gauss-Kronrod quadrature rule to use.
   */
  enum QuadratureRule
  {
    GaussKronrod15 = 1, /**< Use  7-15 points. */
    GaussKronrod21 = 2, /**< Use 10-21 points. */
    GaussKronrod31 = 3, /**< Use 15-31 points. */
    GaussKronrod41 = 4, /**< Use 20-41 points. */
    GaussKronrod51 = 5, /**< Use 25-51 points. */
    GaussKronrod61 = 6  /**< Use 30-61 points. */
  };

  /**
   * \brief Prepares an Integrator for a call to a quadrature function.
   *
   * \param maxSubintervals The maximum number of subintervals allowed in the subdivision process
   *        of quadrature functions. This corresponds to the amount of memory allocated for said
   *        functions.
   */
  Integrator(const int maxSubintervals)
  : m_maxSubintervals(maxSubintervals)
  {
    assert(maxSubintervals >= 1); // \todo use Eigen assert.

    m_errorListIndices.resize(maxSubintervals, 1);
    m_lowerList.resize(maxSubintervals, 1);
    m_upperList.resize(maxSubintervals, 1);
    m_integralList.resize(maxSubintervals, 1);
    m_errorList.resize(maxSubintervals, 1);

    QuadratureKronrod<Scalar>::ComputeNodesAndWeights();
  }

  /**
   * \brief This function calculates an approximation I' to a given definite integral I, the
   *        integral of f from lowerLimit to upperLimit, hopefully satisfying
   *        abs(I - I') <= max(desiredAbsoluteError, desiredRelativeError * abs(I)).
   *
   * This function is best suited for integrands without singularities or discontinuities, which
   * are too difficult for non-adaptive quadrature, and, in particular, for integrands with
   * oscillating behavior of a non-specific type.
   *
   * \param f The function defining the integrand function.
   * \param lowerLimit The lower limit of integration.
   * \param upperLimit The upper limit of integration.
   * \param desiredAbsoluteError The absolute accuracy requested.
   * \param desiredRelativeError The relative accuracy requested. If desiredAbsoluteError <= 0
   *        and desiredRelativeError < max(50*machinePrecision, 0.5e-28), the routine will end with
   *        errorCode = 6.
   * \param quadratureRule The local Gauss-Kronrod quadrature rule to use.
   *
   * \returns The approximation to the integral.
   */
  template <typename FunctionType>
  Scalar quadratureAdaptive(
    const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
    const Scalar desiredAbsoluteError = 1.e-8, const Scalar desiredRelativeError = 1.e-8,
    const QuadratureRule quadratureRule = 1)
  {
    if ((desiredAbsoluteError <= 0.
         && desiredRelativeError < (std::max)(Eigen::NumTraits<Scalar>::epsilon() * Scalar(50.), Scalar(5.e-28) ))
        || m_maxSubintervals < 1)
    {
      m_errorCode = 6;
      return 0.;
    }

    m_errorCode = 0;
    m_numEvaluations = 0;
    m_lowerList[0] = lowerLimit;
    m_upperList[0] = upperLimit;
    m_integralList[0] = 0.;
    m_errorList[0] = 0.;
    m_errorListIndices[0] = 0;
    m_errorListIndices[1] = 1;

    Scalar defAbs;
    Scalar resAbs;

    // First approximation to the integral
    Scalar integral = quadratureKronrod(
      f, lowerLimit, upperLimit, m_estimatedError, defAbs, resAbs, quadratureRule);

    m_numSubintervals = 1;
    m_integralList[0] = integral;
    m_errorList[0] = m_estimatedError;

    // Test on accuracy.
    Scalar errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * Abs(integral));

    if (m_maxSubintervals == 1)
    {
      m_errorCode = 1;
    }
    else if (m_estimatedError <= Eigen::NumTraits<Scalar>::epsilon() * 50. * defAbs
        && m_estimatedError > errorBound)
    {
      m_errorCode = 2;
    }

    if (m_errorCode != 0
        || (m_estimatedError <= errorBound && m_estimatedError != resAbs)
        || m_estimatedError == 0.)
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
    int nrMax = 0;
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
      Scalar defAb1;
      Scalar defAb2;

      const Scalar area1 = quadratureKronrod(
        f, lower1, upper1, error1, resAbs, defAb1, quadratureRule);
      const Scalar area2 = quadratureKronrod(
        f, lower2, upper2, error2, resAbs, defAb2, quadratureRule);

      // Improve previous approximations to integral and error and test for accuracy.
      ++(m_numEvaluations);
      const Scalar area12 = area1 + area2;
      const Scalar error12 = error1 + error2;
      errorSum += error12 - errorMax;
      area += area12 - m_integralList[maxErrorIndex];

      if (defAb1 != error1 && defAb2 != error2)
      {
          if (Abs(m_integralList[maxErrorIndex] - area12) <= Abs(area12) * Scalar(1.e-5)
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

      errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * Abs(area));

      if (errorSum > errorBound)
      {
        // Set error flag in the case that the number of subintervals has reached the max allowable.
        if (m_numSubintervals == m_maxSubintervals)
        {
          m_errorCode = 1;
        }
        // Test for roundoff error and set error flag.
        else if (roundOff1 >= 6 || roundOff2 >= 20)
        {
          m_errorCode = 2;
        }
        // Set m_error_code in the case of poor integrand behaviour within
        // the integration range.
        else if ((std::max)(Abs(lower1), Abs(upper2))
            <= (Eigen::NumTraits<Scalar>::epsilon() * Scalar(100.) + Scalar(1.))
            * (Abs(lower2) + (std::numeric_limits<Scalar>::min)() * Scalar(1.e3) ))
        {
          m_errorCode = 3;
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
      quadratureSort(maxErrorIndex, errorMax,nrMax);

      if (m_errorCode != 0 || errorSum <= errorBound)
      {
        break;
      }
    }

    integral = 0.;

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

  /**
   * \brief Returns the estimated absolute error from the last integration.
   *
   * The value returned will only be valid after calling quadratureAdaptive at least once.
   */
  inline Scalar estimatedError() const {return m_estimatedError;}

  /**
   * \brief Returns the error code.
   *
   * The value returned will only be valid after calling quadratureAdaptive at least once.
   */
  inline int errorCode() const {return m_errorCode;}

private:

  /**
   * \brief This routine maintains the descending ordering in the list of the local error
   *        estimates resulting from the interval subdivision process.
   *
   * At each call two error estimates are inserted using the sequential search method, top-down
   * for the largest error estimate and bottom-up for the smallest error estimate.
   *
   * \param maxErrorIndex The index to the nrMax-th largest error estimate currently in the list.
   * \param errorMax The nrMax-th largest error estimate. errorMaxIndex = errorList(maxError).
   * \param nrMax The integer value such that maxError = errorListIndices(nrMax).
   */
  void quadratureSort(int& maxErrorIndex, Scalar& errorMax, int& nrMax)
  {
    if (m_numSubintervals <= 2)
    {
      m_errorListIndices[0] = 0;
      m_errorListIndices[1] = 1;
      maxErrorIndex = m_errorListIndices[nrMax];
      errorMax = m_errorList[maxErrorIndex];
      return;
    }

    // This part of the routine is only executed if, due to a difficult integrand, subdivision has
    // increased the error estimate. In the normal case the insert procedure should start after the
    // nrMax-th largest error estimate.
    int i = 0;
    int succeed = 0;
    const Scalar errorMaximum = m_errorList[maxErrorIndex];

    if (nrMax != 1)
    {
      for (i = 1; i < nrMax; ++i)
      {
        succeed = m_errorListIndices[nrMax - 1];
        std::cout << "Here 1" << std::endl;

        if (errorMaximum <= m_errorList[succeed])
        {
          std::cout << "Here 2" << std::endl;
          break;
        }

        m_errorListIndices[nrMax] = succeed;
        --nrMax;
      }
    }

    // Compute the number of elements in the list to be maintained in descending order. This number
    // depends on the number of subdivisions remaining allowed.
    int topBegin = m_numSubintervals - 1;
    int bottomEnd = topBegin - 1;
    int start = nrMax + 1;

    if (m_numSubintervals > m_maxSubintervals / 2 + 2)
    {
      topBegin = m_maxSubintervals + 3 - m_numSubintervals + 1;
    }

    // Insert errorMax by traversing the list top-down, starting comparison from the element
    // errorlist(m_errorListIndices(nrMax+1)).
    if (start <= bottomEnd)
    {
      for (i = start; i <= bottomEnd; ++i)
      {
        succeed = m_errorListIndices[i];

        if (errorMaximum >= m_errorList[succeed])
        {
          break;
        }

        m_errorListIndices[i - 1] = succeed;
      }
    }

    if (start > bottomEnd)
    {
      m_errorListIndices[bottomEnd] = maxErrorIndex;
      m_errorListIndices[topBegin] = m_numSubintervals - 1;
      maxErrorIndex = m_errorListIndices[nrMax];
      errorMax = m_errorList[maxErrorIndex];
      return;
    }

    // Insert errorMin by traversing the list bottom-up.
    m_errorListIndices[i - 1] = maxErrorIndex;

    for (int j = i; j <= bottomEnd; ++j)
    {
      succeed = m_errorListIndices[bottomEnd];

      if (m_errorList[m_numSubintervals - 1] < m_errorList[succeed])
      {
        m_errorListIndices[bottomEnd + 1] = m_numSubintervals - 1;
        maxErrorIndex = m_errorListIndices[nrMax];
        errorMax = m_errorList[maxErrorIndex];
        return;
      }

      m_errorListIndices[bottomEnd + 1] = succeed;
      --bottomEnd;
    }

    m_errorListIndices[i] = m_numSubintervals - 1;
    maxErrorIndex = m_errorListIndices[nrMax];
    errorMax = m_errorList[maxErrorIndex];
    return;

  }

  /**
   * \brief This function calculates an approximation I' to a given definite integral I, the
   *        integral of f from lowerLimit to upperLimit and provides an error estimate.
   *
   * \param[in] f The variable representing the function f(x) be integrated.
   * \param[in] lowerLimit The lower limit of integration.
   * \param[in] upperLimit The upper limit of integration.
   * \param[out] errorEstimate Estimate of the modulus of the absolute error, not to exceed
   *             abs(I - I').
   * \param[out] absIntegral The approximation to the integral of abs(f) from lowerLimit to
   *             upperLimit.
   * \param[out] absDiffIntegral The approximation to the integral of
   *             abs(f - I/(upperLimit - lowerLimit)).
   *
   * \returns The approximation I' to the integral I. It is computed by applying the 15, 21, 31,
   *          41, 51, or 61-point kronrod rule obtained by optimal addition of abscissae to the 7,
   *          10, 15, 20, 25, or 30-point Gauss rule.
   *
   * \detail This series of functions represents a priority queue data structure:
   *         - Apply Gauss-Kronrod on the whole initial interval and estimate the error.
   *         - Split the interval symmetrically in two and again apply Gauss-Kronrod on the intervals,
   *           estimate the errors and add a pair (or tuple) "interval, error" to the priority queue.
   *         - If the total error is smaller than the requested tolerance, or if the maximum number
   *           of subdivisions is reached, discontinue the process.
   *         - Otherwise, pop the top element from the queue (highest error), split the interval in two,
   *           and apply Gauss-Kronrod to the new intervals, add those two new elements to the priority
   *           queue and repeat from the previous step.
   */
  template <typename FunctionType>
  Scalar quadratureKronrod(
    const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
    Scalar& estimatedError, Scalar& absIntegral, Scalar& absDiffIntegral,
    const QuadratureRule quadratureRule)
  {
    //QuadratureKronrod<Scalar>::ComputeNodesAndWeights();
    switch (quadratureRule)
    {
    case GaussKronrod15:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod15,
        QuadratureKronrod<Scalar>::weightsGaussKronrod15, QuadratureKronrod<Scalar>::weightsGauss15,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    case GaussKronrod21:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod21,
        QuadratureKronrod<Scalar>::weightsGaussKronrod21, QuadratureKronrod<Scalar>::weightsGauss21,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    case GaussKronrod31:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod31,
        QuadratureKronrod<Scalar>::weightsGaussKronrod31, QuadratureKronrod<Scalar>::weightsGauss31,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    case GaussKronrod41:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod41,
        QuadratureKronrod<Scalar>::weightsGaussKronrod41, QuadratureKronrod<Scalar>::weightsGauss41,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    case GaussKronrod51:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod51,
        QuadratureKronrod<Scalar>::weightsGaussKronrod51, QuadratureKronrod<Scalar>::weightsGauss51,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    case GaussKronrod61:
      return quadratureKronrodHelper(
        QuadratureKronrod<Scalar>::abscissaeGaussKronrod61,
        QuadratureKronrod<Scalar>::weightsGaussKronrod61, QuadratureKronrod<Scalar>::weightsGauss61,
        f, lowerLimit, upperLimit, estimatedError, absIntegral, absDiffIntegral, quadratureRule);

    default:
      return 0.;
    }
  }


  //template <typename FunctionType, int rows1, int rows2>
  template <typename FunctionType, int numKronrodRows, int numGaussRows>
  Scalar quadratureKronrodHelper(
    //some differnces between the template names in IntegratorDefinitions.hpp
    //Array<Scalar, rows1, 1> abscissaeGaussKronrod, Array<Scalar, rows1, 1> weightsGaussKronrod,
    //Array<Scalar, rows2, 1> weightsGauss, const FunctionType& f, const Scalar lowerLimit,
    //const Scalar upperLimit, Scalar& estimatedError, Scalar& absIntegral, Scalar& absDiffIntegral,
    Array<Scalar, numKronrodRows, 1> abscissaeGaussKronrod,
    Array<Scalar, numKronrodRows, 1> weightsGaussKronrod,
    Array<Scalar, numGaussRows, 1> weightsGauss, const FunctionType& f, const Scalar lowerLimit,
    const Scalar upperLimit, Scalar& estimatedError, Scalar& absIntegral, Scalar& absDiffIntegral,
    const QuadratureRule quadratureRule)
  {
    // Half-length of the interval.
    const Scalar halfLength = (upperLimit - lowerLimit) * Scalar(.5);

    // Midpoint of the interval.
    const Scalar center = (lowerLimit + upperLimit) * Scalar(.5);

    const Scalar fCenter = f(center);

    Array<Scalar, numKronrodRows - 1, 1> f1Array;
    Array<Scalar, numKronrodRows - 1, 1> f2Array;

    // The result of the Gauss formula.
    Scalar resultGauss;
    if (quadratureRule % 2 == 0)
    {
      resultGauss = Scalar(0.);
    }
    else
    {
      resultGauss = weightsGauss[weightsGauss.size() - 1] * fCenter;
    }

    // The result of the Kronrod formula.
    Scalar resultKronrod = weightsGaussKronrod[weightsGaussKronrod.size() - 1] * fCenter;
    absIntegral = Abs(resultKronrod);

    for (DenseIndex j = 1; j < weightsGaussKronrod.size() - weightsGauss.size(); ++j)
    {
      const DenseIndex jj = j * 2 - 1;
      const Scalar abscissa = halfLength * abscissaeGaussKronrod[jj];

      const Scalar f1 = f(center - abscissa);
      const Scalar f2 = f(center + abscissa);

      f1Array[jj] = f1;
      f2Array[jj] = f2;
      const Scalar funcSum = f1 + f2;

      resultGauss += weightsGauss[j - 1] * funcSum;
      resultKronrod += weightsGaussKronrod[jj] * funcSum;

      absIntegral += weightsGaussKronrod[jj] * (Abs(f1) + Abs(f2));
    }

    for (DenseIndex j = 0; j < weightsGauss.size(); ++j)
    {
      const DenseIndex jj = j * 2;

      const Scalar abscissa = halfLength * abscissaeGaussKronrod[jj];

      const Scalar f1 = f(center - abscissa);
      const Scalar f2 = f(center + abscissa);

      f1Array[jj] = f1;
      f2Array[jj] = f2;
      const Scalar funcSum = f1 + f2;

      resultKronrod += weightsGaussKronrod[jj] * funcSum;

      absIntegral += weightsGaussKronrod[jj] * (Abs(f1) + Abs(f2));
    }

    // Approximation to the mean value of f over the interval (lowerLimit, upperLimit),
    // i.e. I / (upperLimit - lowerLimit)
    Scalar resultMeanKronrod = resultKronrod * Scalar(.5);

    absDiffIntegral = weightsGaussKronrod[7] * (Abs(fCenter - resultMeanKronrod));

    DenseIndex size1 = weightsGaussKronrod.size() - 1;

    absDiffIntegral += (((f1Array.head(size1) - resultMeanKronrod).abs()
                         + (f2Array.head(size1) - resultMeanKronrod).abs())
                        * weightsGaussKronrod.head(size1)).sum();

    Scalar result = resultKronrod * halfLength;
    absIntegral *= Abs(halfLength);
    absDiffIntegral *= Abs(halfLength);
    estimatedError = Abs((resultKronrod - resultGauss) * halfLength);

    if (absDiffIntegral != Scalar(0.) && estimatedError != Scalar(0.))
    {
      estimatedError = absDiffIntegral
        * (std::min)(Scalar(1.), Pow((estimatedError * Scalar(200.) / absDiffIntegral), Scalar(1.5) ));
    }

    if (absIntegral
        > (std::numeric_limits<Scalar>::min)() / (Eigen::NumTraits<Scalar>::epsilon() * Scalar(50.) ))
    {
      estimatedError = (std::max)(
        Eigen::NumTraits<Scalar>::epsilon() * static_cast<Scalar>(50.) * absIntegral,
        estimatedError);
    }

    return result;
  }


  /**
   * \brief An Array of dimension m_maxSubintervals for error estimates.
   *
   * The first k elements are indices to the error estimates over the subintervals, such that
   * errorList(errorListIndices(0)), ..., errorList(errorListIndices(k - 1)) forms a decreasing
   * sequence, with k = m_numSubintervals if m_numSubintervals <= (m_maxSubintervals/2 + 2),
   * otherwise k = m_maxSubintervals + 1 - m_numSubintervals.
   */
  Array<int, Dynamic, 1> m_errorListIndices;

  /**
   * \brief An Array of dimension m_maxSubintervals for subinterval left endpoints.
   *
   * The first m_numSubintervals elements are the lower end points of the subintervals in the
   * partition of the given integration range (lowerLimit, upperLimit).
   */
  Array<Scalar, Dynamic, 1> m_lowerList;

  /**
   * \brief An Array of dimension m_maxSubintervals for subinterval upper endpoints.
   *
   * The first m_numSubintervals elements are the upper end points of the subintervals in the
   * partition of the given integration range (lowerLimit, upperLimit).
   */
  Array<Scalar, Dynamic, 1> m_upperList;

  /**
   * \brief An Array of dimension m_maxSubintervals for integral approximations.
   *
   * The first m_numSubintervals elements are the integral approximations on the subintervals.
   */
  Array<Scalar, Dynamic, 1> m_integralList;

  /**
   * \brief An Array of dimension m_maxSubintervals for error estimates.
   *
   * The first m_numSubintervals elements of which are the moduli of the absolute error estimates
   * on the subintervals.
   */
  Array<Scalar, Dynamic, 1> m_errorList;

  /**
   * \brief Gives an upper bound on the number of subintervals. Must be at least 1.
   */
  int m_maxSubintervals;

  /**
   * \brief Estimate of the modulus of the absolute error, which should equal or exceed abs(I - I').
   */
  Scalar m_estimatedError;

  /**
   * \brief The number of integrand evaluations.
   */
  int m_numEvaluations;

  /**
   * \brief Error messages generated by the routine.
   *
   * errorCode = 0 Indicates normal and reliable termination of the routine. (It is assumed that
   *               the requested accuracy has been achieved.)
   * errorCode > 0 Indicates abnormal termination of the routine. (The estimates for integral and
   *               m_estimatedError are less reliable and the requested accuracy has not been
   *               achieved.)
   * errorCode = 1 The maximum number of subdivisions allowed has been achieved. One can allow more
   *               subdivisions by increasing the value of m_maxSubintervals. However, if this
   *               yields no improvement it is advised to analyze the integrand in order to
   *               determine the integration difficulaties. If the position of a local difficulty
   *               can be determined, (i.e. singularity or discontinuity within the interval), one
   *               will probably gain from splitting up the interval at this point and calling the
   *               integrator on the subranges. If possible, an appropriate special-purpose
   *               integrator should be used which is designed for handling the type of difficulty
   *               involved.
   * errorCode = 2 The occurrence of roundoff error is detected, preventing the requested tolerance
   *               from being achieved.
   * errorCode = 3 Extremely bad integrand behaviour occurs at points in the integration interval.
   * errorCode = 6 The input is invalid, because (desiredAbsoluteError <= 0 and
   *               desiredRealtiveError < 50 * relativeMachineAccuracy, or
   *               m_maxSubintervals < 1. \todo make relativeMachineAccuracy a member variable.
   */
  int m_errorCode;

  /**
   * \brief The number of subintervals actually produced in the subdivision process.
   */
   int m_numSubintervals;

};

}

//#include "IntegratorDefinitions.hpp"

#endif // EIGEN_INTEGRATOR_H
