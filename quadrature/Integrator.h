/**
 * \file Integrator.h
 * \sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package
 *     for Automatic Integration, Springer Verlag, 1983.
 */

#ifndef EIGEN_INTEGRATOR_H
#define EIGEN_INTEGRATOR_H

#include <Eigen/Dense>

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
  Integrator(const int maxSubintervals);

  /**
   * \brief This function calculates an approximation I' to a given definite integral I, the
   *        integral of f from lowerLimit to upperLimit, hopefully satisfying
   *        fabs(I - I') <= max(desiredAbsoluteError, desiredRelativeError * abs(I)).
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
    const QuadratureRule quadratureRule = 1);

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
  void quadratureSort(int& maxErrorIndex, Scalar& errorMax, int& nrMax);

  /**
   * \brief This function calculates an approximation I' to a given definite integral I, the
   *        integral of f from lowerLimit to upperLimit and provides an error estimate.
   *
   * \param[in] f The variable representing the function f(x) be integrated.
   * \param[in] lowerLimit The lower limit of integration.
   * \param[in] upperLimit The upper limit of integration.
   * \param[out] errorEstimate Estimate of the modulus of the absolute error, not to exceed
   *             fabs(I - I').
   * \param[out] absIntegral The approximation to the integral of fabs(f) from lowerLimit to
   *             upperLimit.
   * \param[out] absDiffIntegral The approximation to the integral of
   *             fabs(f - I/(upperLimit - lowerLimit)).
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
    const QuadratureRule quadratureRule);

  template <typename FunctionType, int rows1, int rows2>
  Scalar quadratureKronrodHelper(
    Array<Scalar, rows1, 1> abscissaeGaussKronrod, Array<Scalar, rows1, 1> weightsGaussKronrod,
    Array<Scalar, rows2, 1> weightsGauss, const FunctionType& f, const Scalar lowerLimit,
    const Scalar upperLimit, Scalar& estimatedError, Scalar& absIntegral, Scalar& absDiffIntegral,
    const QuadratureRule quadratureRule);

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
   * \brief Estimate of the modulus of the absolute error, which should equal or exceed fabs(I - I').
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
   *               desiredRealtiveError < max(50*relativeMachineAccuracy, 0.5e-28)), or
   *               m_maxSubintervals < 1. \todo make relativeMachineAccuracy a member variable.
   */
  int m_errorCode;

  /**
   * \brief The number of subintervals actually produced in the subdivision process.
   */
   int m_numSubintervals;

};

}

#include "IntegratorDefinitions.h"

#endif // EIGEN_INTEGRATOR_H
