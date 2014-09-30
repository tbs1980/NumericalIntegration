/**
 * \file Integrator.cpp
 * \sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package
 *     for Automatic Integration, Springer Verlag, 1983.
 */

#ifndef EIGEN_INTEGRATOR_DEFINITIONS_H
#define EIGEN_INTEGRATOR_DEFINITIONS_H

//#include "QuadratureKronrod.h"
#include "QuadratureKronrod_STB.hpp"
#include <iostream>

namespace Eigen
{

/**
 * \todo This is inlined so it can reside in header. It won't be when the class is templated for
 *       arbitrary scalar types.
 */
template <typename Scalar_>
inline Integrator<Scalar_>::Integrator(const int maxSubintervals) : m_maxSubintervals(maxSubintervals)
{
  assert(maxSubintervals >= 1); // \todo use Eigen assert.

  m_errorListIndices.resize(maxSubintervals, 1);
  m_lowerList.resize(maxSubintervals, 1);
  m_upperList.resize(maxSubintervals, 1);
  m_integralList.resize(maxSubintervals, 1);
  m_errorList.resize(maxSubintervals, 1);
}

/**
 * \todo This is inlined so it can reside in header. It won't be when the class is templated for
 *       arbitrary scalar types.
 */
template <typename Scalar_>
inline void Integrator<Scalar_>::quadratureSort(int& maxErrorIndex, Scalar& errorMax, int& nrMax)
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

template <typename Scalar_>
template <typename FunctionType>
Scalar_ Integrator<Scalar_>::quadratureAdaptive(
  const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
  const Scalar desiredAbsoluteError, const Scalar desiredRelativeError,
  const QuadratureRule quadratureRule)
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
  Scalar errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * std::abs(integral));

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
        if (std::abs(m_integralList[maxErrorIndex] - area12) <= std::abs(area12) * Scalar(1.e-5)
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

    errorBound = (std::max)(desiredAbsoluteError, desiredRelativeError * std::abs(area));

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
      else if ((std::max)(std::abs(lower1), std::abs(upper2))
          <= (Eigen::NumTraits<Scalar>::epsilon() * Scalar(100.) + Scalar(1.))
          * (std::abs(lower2) + (std::numeric_limits<Scalar>::min)() * Scalar(1.e3) ))
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
    std::cout << "maxErrorIndex" << maxErrorIndex << "  :  ";
    quadratureSort(maxErrorIndex, errorMax,nrMax);
    std::cout << maxErrorIndex << std::endl;
    //errorMax = m_errorList.maxCoeff(&maxErrorIndex);

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

template <typename Scalar_>
template <typename FunctionType>
Scalar_ Integrator<Scalar_>::quadratureKronrod(
  const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
  Scalar& estimatedError, Scalar& absIntegral,  Scalar& absDiffIntegral,
  const QuadratureRule quadratureRule)
{
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

template <typename Scalar_>
template <typename FunctionType, int numKronrodRows, int numGaussRows>
Scalar_ Integrator<Scalar_>::quadratureKronrodHelper(
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
  absIntegral = std::abs(resultKronrod);

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

    absIntegral += weightsGaussKronrod[jj] * (std::abs(f1) + std::abs(f2));
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

    absIntegral += weightsGaussKronrod[jj] * (std::abs(f1) + std::abs(f2));
  }

  // Approximation to the mean value of f over the interval (lowerLimit, upperLimit),
  // i.e. I / (upperLimit - lowerLimit)
  Scalar resultMeanKronrod = resultKronrod * Scalar(.5);

  absDiffIntegral = weightsGaussKronrod[7] * (std::abs(fCenter - resultMeanKronrod));

  DenseIndex size1 = weightsGaussKronrod.size() - 1;

  absDiffIntegral += (((f1Array.head(size1) - resultMeanKronrod).abs()
                       + (f2Array.head(size1) - resultMeanKronrod).abs())
                      * weightsGaussKronrod.head(size1)).sum();

  Scalar result = resultKronrod * halfLength;
  absIntegral *= std::abs(halfLength);
  absDiffIntegral *= std::abs(halfLength);
  estimatedError = std::abs((resultKronrod - resultGauss) * halfLength);

  if (absDiffIntegral != Scalar(0.) && estimatedError != Scalar(0.))
  {
    estimatedError = absDiffIntegral
      * (std::min)(Scalar(1.), std::pow((estimatedError * Scalar(200.) / absDiffIntegral), Scalar(1.5) ));
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

}

#endif // EIGEN_INTEGRATOR_DEFINITIONS_H
