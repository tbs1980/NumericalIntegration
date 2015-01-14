/**
 * \file
 * \brief - The routine calculates an approximation integral to a given Fourier integral
 *           i = integeral of f(x)*w(x) over (lowerLimit,infiniteBoundKey) where w(x) = cos(omega*x) or w(x) = sin(omega*x),
 *           hopefully satisfying following claim for accuracy abs(i-integral) <= desiredAbsoluteError.
 *
 * \details This routine is a special-purpose automatic integerator for integeration between zeros with dqawoe with convergence acceleration with dqelg. 
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] omega  - parameter in the weight function
 * \param[] integer - Indicates which weight function is to be used: if integer = 1, w(x) = cos(omega*x), and  if integer = 2, w(x) = sin(omega*x)
 *             If integer != 1 && integer != 2, the routine will end with m_errorCode = 6.
 * \param[] desiredAbsoluteError - Absolute accuracy requested, desiredAbsoluteError > 0.  If desiredAbsoluteError <= 0, the routine will end with m_errorCode = 6.
 * \param[] limlist - limlist gives an upper finiteBound on the number of cycles, limlist >= 1. If limlist < 3, the routine will end with m_errorCode = 6.
 * \param[] m_maxSubintervals - Gives an upper finiteBound on the number of subintervals allowed in the partition of each cycle, m_maxSubintervals >= 1 each cycle, m_maxSubintervals >= 1.
 * \param[] maxp1 - gives an upper finiteBound on the number of Chebyshev moments which can be stored, i.e. for the intervals of lengths abs(upperLimit-lowerLimit)*2**(-l), l=0,1, ..., maxp1-2, maxp1 >= 1
 * \param[] rslist - Vector of dimension at least limlist rslist(k) contains the integeral contribution over the interval (a+(k-1)c,a+kc) where c = (2*int(abs(omega))+1)*M_PI/abs(omega), from k = 1, to list.  Note that, if omega = 0, rslist(1) contains the value of the integeral over (a,infiniteBoundKey).
 * \param[] erlist - vector of dimension at least limlist erlist(k) contains the error estimate corresponding with rslist(k).
 * \param[] list - number of subintervals needed for the integeration.  If omega = 0 then list is set to 1.
 * \param[] chebyshevMoments - array of dimension at least (maxp1,25), providing space for the chebyshev moments needed within the cycles
 * \param[] c1, c2 - end points of subinterval (of length cycle)
 * \param[] cycle - (2*int(abs(omega))+1)*M_PI/abs(omega)
 * \param[] partialIntegral - vector of dimension at least (limexp+2) (see routine qelg) partialIntegral contains the part of the epsilon table which is still needed for further computations. Each element of partialIntegral is a partial sum of the series which should sum to the value of the integeral.
 * \param[] errorSum - sum of error estimates over the subintervals, calculated cumulatively
 * \param[] desiredAbsoluteError - absolute tolerance requested over current subinterval
 * \param[] chebyshevMoments - array containing the modified chebyshev moments (see also routine qc25f)
 *
 * \returns The approximation to the integeral.
 */

template <typename FunctionType>
Scalar adaptiveQuadratureForFourierIntegrals(
        const FunctionType& f, const Scalar lowerLimit, const Scalar upperLimit,
        const Scalar desireabsoluteError = Scalar(0.), const Scalar desiredRelativeError = Scalar(0.),
        const QuadratureRule quadratureRule = 1)
{
    if ((desiredAbsoluteError <= Scalar(0.) && desiredRelativeError < Eigen::NumTraits<Scalar>::epsilon())
        || (integer != 1 && integer != 2)
        || limlist < 3 
        || leniw < (limlist+2) 
        || lenw < (leniw * 2 + maxp1 * 25)
        || m_maxSubintervals < 1
        || maxp1 < 1 
        || m_errorCode !0)

    {
        std::cout << "Abnormal return.  m_errorCode = " << m_errorCode << std::endl;
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
    Scalar lowerLimit;
    Scalar m_estimatedError;
    Scalar desiredAbsoluteError;
    Scalar omega;
    Scalar integral;
    Scalar work;

    int m_errorCode;
    int integer;
    int leniw;
    int limlist;
    int list;
    int maxp1;

    iwork(leniw);
    work(lenw);

    Scalar integral = 0.;
    
    m_errorCode = 0;
    m_maxSubintervals = (leniw-limlist) / 2;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    m_estimatedError = 0.;

    int l1 = limlist + 1;
    int l2 = limlist + l1;
    int l3 = m_maxSubintervals+l2;
    int l4 = m_maxSubintervals+l3;
    int l5 = m_maxSubintervals+l4;
    int l6 = m_maxSubintervals+l5;
    int ll2 = m_maxSubintervals+l1;
    int lvl = 0;

    real a,absoluteEpsilon,m_estimatedError,m_lowerList,m_upperList,chebyshevMoments,errorCorrection,cycle,c1,c2,dl,dla,drl,list,ep,eps,desiredAbsoluteError,desiredAbsoluteError,erlist,errorSum,fact,omega,p,M_PI,p1,partialIntegral,resultEpsilon,integral,res3la,m_integralList,rslist,r1mach,(std::numeric_limits<Scalar>::min)()

    int ktmin,
    int l,
    int list = 0;
    ll,
    maxp1,
    nEval,
    nnlog,
    numberOfExtrapolations,
    numberElementsList2

    chebyshevMoments(maxp1,25);
    erlist(limlist);
    nnlog(m_maxSubintervals);
    partialIntegral(52);
    res3la(3);
    rslist(limlist);

    // The dimension of  partialIntegral  is determined by the value of limexp in subroutine qelg (partialIntegral must be of dimension (limexp+2) at least).
    Scalar p / 0.9;

    integral = 0.;
    m_estimatedError = 0.;
    m_numEvaluations = 0;
    m_errorCode = 0;

    if(omega == 0.)
    {
        // integeration by qagie if omega is zero
        if(integer == 1) 
        {
            call qagie(f,0.,1,desiredAbsoluteError,0.,m_maxSubintervals,integral,m_estimatedError,m_numEvaluations,m_errorCode,m_lowerList,m_upperList,m_integralList,list,listIndices,m_numSubintervals)
        }

        rslist(1) = integral;
        erlist(1) = m_estimatedError;
        m_errorCodeList(1) = m_errorCode;
        list = 1;
        return integral;
    }

   10 l = abs(omega)
      dl = 2*l+1
      cycle = dl*M_PI/abs(omega)
      m_errorCode = 0
      ktmin = 0
      m_numEvaluations = 0
      numberElementsList2 = 0
      numberOfExtrapolations = 0
      c1 = a
      c2 = cycle+a
      p1 = 1.-p
      eps = desiredAbsoluteError
      (std::numeric_limits<Scalar>::min)() = r1mach(1)
      if(desiredAbsoluteError > (std::numeric_limits<Scalar>::min)()/p1) eps = desiredAbsoluteError*p1
      ep = eps
      fact = 1.
      errorCorrection = 0.
      m_estimatedError = 0.
      errorSum = 0.

  // Main loop for the integration
      do 50 list = 1,limlist

  // integerate over current subinterval.
        dla = list
        desiredAbsoluteError = eps*fact
        call qawoe(f,c1,c2,omega,integer,desiredAbsoluteError,0.,m_maxSubintervals,list,maxp1,rslist(list),erlist(list),nEval,m_errorCodeList(list),m_numSubintervals,m_lowerList,m_upperList,m_integralList,list,listIndices,nnlog,momentsComputed,chebyshevMoments)
        m_numEvaluations = m_numEvaluations+nEval
        fact = fact*p
        errorSum = errorSum+erlist(list)
        drl = 5.*abs(rslist(list))

    // test on accuracy with partial sum
        if(errorSum+drl <= desiredAbsoluteError && list >= 6) go to 80
        errorCorrection = (std::max)(errorCorrection,erlist(list))
        if(m_errorCodeList(list) != 0) eps = (std::max)(ep,errorCorrection*p1)
        if(m_errorCodeList(list) != 0) m_errorCode = 7
        if(m_errorCode == 7 && (errorSum+drl) <= errorCorrection*10. && list > 5) go to 80
        numberElementsList2 = numberElementsList2+1
        if(list > 1) go to 20
        partialIntegral(1) = rslist(1)
        go to 40
   20   partialIntegral(numberElementsList2) = partialIntegral(ll)+rslist(list)
        if(list == 2) go to 40

    // test on maximum number of subintervals
        if(list == limlist) m_errorCode = 1

    // perform new extrapolationPerformed
        call qelg(numberElementsList2,partialIntegral,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)

    // test whether extrapolationPerformedolated integral is infiniteBoundKeyluenced by roundoff
        ktmin = ktmin+1
        if(ktmin >= 15 && m_estimatedError <= 0.001*(errorSum+drl)) m_errorCode = 4
        if(absoluteEpsilon > m_estimatedError && list != 3) go to 30
        m_estimatedError = absoluteEpsilon
        integral = resultEpsilon
        ktmin = 0

    // if m_errorCode is not 0, check whether direct integral (partial sum) or extrapolationPerformedolated integral yields the best integeral approximation
        if((m_estimatedError+10.*errorCorrection) <= desiredAbsoluteError || (m_estimatedError <= desiredAbsoluteError && 10.*errorCorrection >= desiredAbsoluteError)) go to 60
   30   if(m_errorCode != 0 && m_errorCode != 7) go to 60
   40   ll = numberElementsList2
        c1 = c2
        c2 = c2+cycle
   50 continue

    // set final integral and error estimate
   60 m_estimatedError = m_estimatedError+10.*errorCorrection
      if(m_errorCode == 0) go to 999
      if(integral != 0. && partialIntegral(numberElementsList2) != 0.) go to 70
      if(m_estimatedError > errorSum) go to 80
      if(partialIntegral(numberElementsList2) == 0.) go to 999
   70 if(m_estimatedError/abs(integral) > (errorSum+drl)/abs(partialIntegral(numberElementsList2))) go to 80
      if(m_errorCode >= 1 && m_errorCode != 7) m_estimatedError = m_estimatedError+drl
      go to 999
   80 integral = partialIntegral(numberElementsList2)
      m_estimatedError = errorSum+drl
  999 return
      end
