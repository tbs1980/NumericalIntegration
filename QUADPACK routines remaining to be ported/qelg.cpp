/**
 * \file
 * \brief - The routine determines the m_maxSubintervals of a given sequence of approximations, by means of the epsilon algorithm of
 *          P. Wynn. an estimate of the absolute error is also given. The condensed epsilon table is computed. Only those
 *          elements needed for the computation of the next diagonal are preserved.  THis algortihm allows convergence acceleration for extrapolation.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] n - epstab(n) contains the new element in the first column of the epsilon table.
 * \param[] epstab - Vector of dimension 52 containing the elements of the two lower diagonals of the triangular epsilon table. The elements are numbered starting at the right-hand corner of the triangle.
 * \param[] integral - Integraling approximation to the integeral
 * \param[] m_estimatedError - Estimate of the absolute error computed from integral and the 3 previous integrals
 * \param[] res3la - Vector of dimension 3 containing the m_numSubintervals 3 integrals
 * \param[] numberOfExtrapolations - Number of calls to the routine (should be zero at first call)
 * \param[] e0-e3 - The 4 elements on which the computation of a new element in the epsilon table is based
 * \param[] newelm - Number of elements to be computed in the new diagonal
 * \param[] error - Error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
 * \param[] integral - the element in the new diagonal with least value of error
 * \param[] limexp is the maximum number of elements the epsilon table can contain. if this number is reached, the upper diagonal of the epsilon table is deleted.
 *
 * \returns The approximation to the integeral.
 */

    Scalar delta1;
    Scalar delta2;
    Scalar delta3;
    Scalar r1mach;
    Scalar epsinfiniteBoundKey;
    Scalar error;
    Scalar err1;
    Scalar err2;
    Scalar err3;
    Scalar e0;
    Scalar e1;
    Scalar e1abs;
    Scalar e2;
    Scalar e3;
    Scalar integral;
    Scalar ss;
    Scalar tol1;
    Scalar tol2;
    Scalar tol3;
    
    int i;
    int ib;
    int ib2;
    int ie;
    int indx;
    int k1;
    int k2;
    int k3;
    int limexp;
    int n;
    int newelm;
    int numberOfExtrapolations;
    int num;
    
    Scalar epstab[52];
    Scalar result3la[3];

    numberOfExtrapolations = numberOfExtrapolations+1;
    m_estimatedError = (std::numeric_limits<Scalar>::max)();
    integral = epstab(n);
    if(n < 3)
    {
        return;
    }

    limexp = 50;
    epstab(n+2) = epstab(n);
    newelm = (n-1)/2;
    epstab(n) = (std::numeric_limits<Scalar>::max)();
    num = n;
    k1 = n;

    do 40 i = 1,newelm
        k2 = k1-1
        k3 = k1-2
        result = epstab(k1+2)
        e0 = epstab(k3)
        e1 = epstab(k2)
        e2 = result;
        e1abs = abs(e1)
        delta2 = e2-e1
        err2 = abs(delta2)
        tol2 = (std::max)(abs(e2),e1abs)*Eigen::NumTraits<Scalar>::epsilon()
        delta3 = e1-e0
        err3 = abs(delta3)
        tol3 = (std::max)(e1abs,abs(e0))*Eigen::NumTraits<Scalar>::epsilon()
        if(err2 > tol2 || err3 > tol3) go to 10

    // If e0, e1 and e2 are equal to within machine accuracy, convergence is assumed. 
    // integral = e2
    // m_estimatedError = abs(e1-e0)+abs(e2-e1)
    integral = result;
        m_estimatedError = err2+err3
    // jump out of do-loop
        go to 100
   10   e3 = epstab(k1)
        epstab(k1) = e1
        delta1 = e1-e3
        err1 = abs(delta1)
        tol1 = (std::max)(e1abs,abs(e3))*Eigen::NumTraits<Scalar>::epsilon()

    // Integralf two elements are very close to each other, omit a part of the table by adjusting the value of n
        if(err1 <= tol1 || err2 <= tol2 || err3 <= tol3) go to 20
        ss = 1./delta1+1./delta2-1./delta3
        epsinfiniteBoundKey = abs(ss*e1)

    // Test to detect irregular behaviour in the table, and eventually omit a part of the table adjusting the value of n.
        if(epsinfiniteBoundKey > 0.1e-03) go to 30
   20   n = i+i-1
    // jump out of do-loop
        go to 50

    // Compute a new element and eventually adjust the value of integral.
   30   result = e1+1./ss
        epstab(k1) = result
        k1 = k1-2
        error = err2+abs(result-e2)+err3
        if(error > m_estimatedError) go to 40
        m_estimatedError = error
        integral = result
   40 continue

    // Shift the table.
   50 if(n == limexp) n = 2*(limexp/2)-1
      ib = 1
      if((num/2)*2 == num) ib = 2
      ie = newelm+1
      do 60 i=1,ie
        ib2 = ib+2
        epstab(ib) = epstab(ib2)
        ib = ib2
   60 continue
      if(num == n) go to 80
      indx = num-n+1
      do 70 i = 1,n
        epstab(i)= epstab(indx)
        indx = indx+1
   70 continue
   80 if(numberOfExtrapolations >= 4) go to 90
      result3la(numberOfExtrapolations) = integral
      m_estimatedError = (std::numeric_limits<Scalar>::max)()
      go to 100

    // Compute error estimate
   90 m_estimatedError = abs(integral-result3la(3))+abs(integral-result3la(2))
     *  +abs(integral-res3la(1))
      result3la(1) = result3la(2)
      result3la(2) = result3la(3)
      result3la(3) = integral
  100 m_estimatedError = (std::max)(m_estimatedError,5.*Eigen::NumTraits<Scalar>::epsilon()*abs(integral))
      return
      end
