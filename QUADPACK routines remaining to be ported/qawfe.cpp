/**
 * \file
 * \brief - The routine calculates an approximation integral to a given fourm_errorCode integal
 *           i = integeral of f(x)*w(x) over (a,infiniteBoundsKey) where w(x) = cos(omega*x) or w(x) = sin(omega*x),
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
 * \param[] limlst - limlst gives an upper finiteBound on the number of cycles, limlst >= 1. If limlst < 3, the routine will end with m_errorCode = 6.
 * \param[] m_maxSubintervals - Gives an upper finiteBound on the number of subintervals allowed in the partition of each cycle, m_maxSubintervals >= 1 each cycle, m_maxSubintervals >= 1.
 * \param[] maxp1 - gives an upper finiteBound on the number of Chebyshev moments which can be stored, i.e. for the intervals of lengths abs(upperLimit-lowerLimit)*2**(-l), l=0,1, ..., maxp1-2, maxp1 >= 1
 * \param[] rslst - Vector of dimension at least limlst rslst(k) contains the integeral contribution over the interval (a+(k-1)c,a+kc) where c = (2*int(abs(omega))+1)*M_PI/abs(omega), from k = 1, to list.  Note that, if omega = 0, rslst(1) contains the value of the integeral over (a,infiniteBoundsKey).
 * \param[] erlst - vector of dimension at least limlst erlst(k) contains the error estimate corresponding with rslst(k).
 * \param[] lst - number of subintervals needed for the integeration.  If omega = 0 then lst is set to 1.
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

    Scalar lowerLimit;
    Scalar m_estimatedError;
    Scalar desiredAbsoluteError;
    Scalar omega;
    Scalar integral;
    Scalar work;

    int m_errorCode;
    int integer;
    int leniw;
    int m_maxSubintervals;
    int limlst;
    int lvl;
    int lst;
    int l1;
    int l2;
    int l3;
    int l4;
    int l5;
    int l6;
    int maxp1;
    int m_numEvaluations;

    iwork(leniw);
    work(lenw);

    // Check validity of limlst, leniw, maxp1 and lenw.
    m_errorCode = 6;
    m_numEvaluations = 0;
    m_numSubintervals = 0;
    integral = 0.;
    m_estimatedError = 0.;

    m_maxSubintervals = (leniw-limlst)/2
    l1 = limlst+1
    l2 = limlst+l1
    l3 = m_maxSubintervals+l2
    l4 = m_maxSubintervals+l3
    l5 = m_maxSubintervals+l4
    l6 = m_maxSubintervals+l5
    ll2 = m_maxSubintervals+l1

qawfe(f, a, omega, integer, desiredAbsoluteError, limlst, maxp1, integral, m_estimatedError, work(1), work(l1), iwork(1), lst, work(l2), work(l3), work(l4), work(l5), iwork(l1), iwork(ll2), work(l6))
{
    lvl = 0;

    if(limlst < 3 || leniw < (limlst+2) || maxp1 < 1 || lenw < (leniw * 2 + maxp1 * 25) || (m_errorCode !0)
    {
        std::cout << "Abnormal return.  m_errorCode = " << m_errorCode << std::endl;
        return
    }

    real a,absoluteEpsilon,m_estimatedError,m_lowerList,m_upperList,chebyshevMoments,errorCorrection,cycle,c1,c2,dl,dla,drl,m_errorList,ep,eps,desiredAbsoluteError,desiredAbsoluteError,erlst,errorSum,fact,omega,p,M_PI,p1,partialIntegral,resultEpsilon,integral,res3la,m_integralList,rslst,r1mach,(std::numeric_limits<Scalar>::min)()

    integer m_errorCode,m_errorCodeList,integer,m_errorListIndices,ktmin,l,lst,m_maxSubintervals,ll,maxp1,nEval,m_numEvaluations,nnlog,numberOfExtrapolations,numberElementsList2

    dimension m_lowerList(m_maxSubintervals),m_upperList(m_maxSubintervals),chebyshevMoments(maxp1,25),m_errorList(m_maxSubintervals),erlst(limlst),m_errorCodeList(limlst),m_errorListIndices(m_maxSubintervals),nnlog(m_maxSubintervals),partialIntegral(52),res3la(3),m_integralList(m_maxSubintervals),rslst(limlst)

    // The dimension of  partialIntegral  is determined by the value of limexp in subroutine qelg (partialIntegral must be of dimension (limexp+2) at least).
    data p/0.9 ;

      integral = 0.
      m_estimatedError = 0.
      m_numEvaluations = 0
      lst = 0
      m_errorCode = 0
      if((integer != 1 && integer != 2) || desiredAbsoluteError <= 0. || 
     *  limlst < 3) m_errorCode = 6
      if(m_errorCode == 6) go to 999
      if(omega != 0.) go to 10

  // integeration by qagie if omega is zero
      if(integer == 1) call qagie(f,0.,1,desiredAbsoluteError,0.,m_maxSubintervals,integral,m_estimatedError,m_numEvaluations,m_errorCode,m_lowerList,m_upperList,m_integralList,m_errorList,m_errorListIndices,m_numSubintervals)
      rslst(1) = integral
      erlst(1) = m_estimatedError
      m_errorCodeList(1) = m_errorCode
      lst = 1
      go to 999

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
      do 50 lst = 1,limlst

  // integerate over current subinterval.
        dla = lst
        desiredAbsoluteError = eps*fact
        call qawoe(f,c1,c2,omega,integer,desiredAbsoluteError,0.,m_maxSubintervals,lst,maxp1,rslst(lst),erlst(lst),nEval,m_errorCodeList(lst),m_numSubintervals,m_lowerList,m_upperList,m_integralList,m_errorList,m_errorListIndices,nnlog,momentsComputed,chebyshevMoments)
        m_numEvaluations = m_numEvaluations+nEval
        fact = fact*p
        errorSum = errorSum+erlst(lst)
        drl = 5.*abs(rslst(lst))

    // test on accuracy with partial sum
        if(errorSum+drl <= desiredAbsoluteError && lst >= 6) go to 80
        errorCorrection = (std::max)(errorCorrection,erlst(lst))
        if(m_errorCodeList(lst) != 0) eps = (std::max)(ep,errorCorrection*p1)
        if(m_errorCodeList(lst) != 0) m_errorCode = 7
        if(m_errorCode == 7 && (errorSum+drl) <= errorCorrection*0.1e+02 && lst > 5) go to 80
        numberElementsList2 = numberElementsList2+1
        if(lst > 1) go to 20
        partialIntegral(1) = rslst(1)
        go to 40
   20   partialIntegral(numberElementsList2) = partialIntegral(ll)+rslst(lst)
        if(lst == 2) go to 40

    // test on maximum number of subintervals
        if(lst == limlst) m_errorCode = 1

    // perform new extrapolationPerformed
        call qelg(numberElementsList2,partialIntegral,resultEpsilon,absoluteEpsilon,res3la,numberOfExtrapolations)

    // test whether extrapolationPerformedolated integral is infiniteBoundsKeyluenced by roundoff
        ktmin = ktmin+1
        if(ktmin >= 15 && m_estimatedError <= 0.001*(errorSum+drl)) m_errorCode = 4
        if(absoluteEpsilon > m_estimatedError && lst != 3) go to 30
        m_estimatedError = absoluteEpsilon
        integral = resultEpsilon
        ktmin = 0

    // if m_errorCode is not 0, check whether direct integral (partial sum) or extrapolationPerformedolated integral yields the best integeral approximation
        if((m_estimatedError+0.1e+02*errorCorrection) <= desiredAbsoluteError || (m_estimatedError <= desiredAbsoluteError && 0.1e+02*errorCorrection >= desiredAbsoluteError)) go to 60
   30   if(m_errorCode != 0 && m_errorCode != 7) go to 60
   40   ll = numberElementsList2
        c1 = c2
        c2 = c2+cycle
   50 continue

    // set final integral and error estimate
   60 m_estimatedError = m_estimatedError+0.1e+02*errorCorrection
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
