/**
 * \file
 * \brief - To compute i = integeral of f*w over (bl,br), with error estimate, where the weight function w has a singular
 *          behaviour of algebraico-logarithmic type at the points lowerLimit and/or upperLimit. (bl,br) is a part of (lowerLimit, upperLimit).
 *          The routine employs 25-point Clenshaw-Curtis Integeration. 
 * \details - This routine contains the integeration rules for integerands having algebraico-logarithmic end point singularities
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] bl - lower m_maxSubintervals of integeration, bl >= a
 * \param[] br - upper m_maxSubintervals of integeration, br <= b
 * \param[] alpha - Parameter in the weight function
 * \param[] beta - Parameter in the weight function
 * \param[] ri,rj,rg,rh - Modified chebyshev moments for the application of the generalized clenshaw-curtis method (computed in subroutine dqmomo)
 * \param[] integral - Approximation to the integeral integral is computed by using a generalized
 *              Clenshaw-Curtis method if b1 = a or br = b. In all other cases the 15-point kronrod
 *              rule is applied, obtained by optimal addition of abscissae to the 7-point gauss rule.
 * \param[] m_estimatedError - estimate of the modulus of the absolute error, which should equal or exceed abs(i-integral)*
 * \param[] m_estimatedError - approximation to the integeral of abs(f*w-i/(upperLimit-lowerLimit))*
 * \param[] integer - which determines the weight function
 *                    = 1   w(x) = (x-a)**alpha*(b-x)**beta
 *                    = 2   w(x) = (x-a)**alpha*(b-x)**beta*log(x-a)
 *                    = 3   w(x) = (x-a)**alpha*(b-x)**beta*log(b-x)
 *                    = 4   w(x) = (x-a)**alpha*(b-x)**beta*log(x-a)*log(b-x)
 * \param[] nEval - Number of integerand evaluations
 * \param[] fValue - value of the function f at the points (br-bl)*0.5*cos(k*M_PI/24)+(br+bl)*0.5 k = 0, ..., 24
 * \param[] chebyshevDegree12 - coefficients of the chebyshev series expansion of degree 12, for the function f, in the interval (bl,br)
 * \param[] chebyshevDegree24 - coefficients of the chebyshev series expansion of degree 24, for the function f, in the interval (bl,br)
 * \param[] integeral12 - approximation to the integeral obtained from chebyshevDegree12
 * \param[] integeral24 - approximation to the integeral obtained from chebyshevDegree24
 * \param[] qwgts - external function subprogram defining the four possible weight functions
 *
 * \returns The approximation to the integeral.
 */

    Scalar lowerLimit,m_estimatedError,alphlowerLimit, upperLimit,betlowerLimit, upperLimitl,br,center,chebyshevDegree12,chebyshevDegree24,dc,f,factor,fix,fValue,halfLength,absIntegral,m_estimatedError,integral,integeral12,integeral24,rg,rh,ri,rj,u,qwgts,x
    int i;
    int integer;
    int isym;
    int nEval

    chebyshevDegree12(13),chebyshevDegree24(25),fValue(25),rg(25),rh(25),ri(25),rj(25),x(11)

    // The vector x contains the values cos(k*M_PI/24) k = 1, ..., 11, to be used for the computation of the Chebyshev series expansion of f.
    data x
    {
        0.9914448613738104,
        0.9659258262890683,
        0.9238795325112868,
        0.8660254037844386,
        0.7933533402912352,
        0.7071067811865475,
        0.6087614290087206,
        0.5000000000000000,
        0.3826834323650898,
        0.2588190451025208,
        0.1305261922200516
    }.finished();

    nEval = 25
    if(bl == lowerLimit && (alpha != 0. || integer == 2 || integer == 4)) go to 10
    if(br == upperLimit && (beta != 0. || integer == 3 || integer == 4)) go to 140

    // If lowerLimit > bl and upperLimit < br, apply the 15-point gauss-kronrod scheme.
    call qk15w(f,qwgts,lowerLimit, upperLimit,alphlowerLimit, beta,integer,bl,br,integral,m_estimatedError,absIntegral,m_estimatedError)
    nEval = 15
    go to 270

    // This part of the program is executed only if a = bl.
    // Compute the chebyshev series expansion of the following function
    // f1 = (0.5*(upperLimit+upperLimit-br-lowerLimit)-0.5*(br-lowerLimit)*x)**beta*f(0.5*(br-lowerLimit)*x+0.5*(br+lowerLimit))
   10 halfLength = 0.5*(br-bl)
      center = 0.5*(br+bl)
      fix = upperLimit-center
      fValue(1) = 0.5*f(halfLength+center)*(fix-halfLength)**beta
      fValue(13) = f(center)*(fix**beta)
      fValue(25) = 0.5*f(center-halfLength)*(fix+halfLength)**beta
      do 20 i=2,12
        u = halfLength*x(i-1)
        isym = 26-i
        fValue(i) = f(u+center)*(fix-u)**beta
        fValue(isym) = f(center-u)*(fix+u)**beta
   20 continue
      factor = halfLength**(alpha+1.)
      integral = 0.
      m_estimatedError = 0.
      integeral12 = 0.
      integeral24 = 0.
      if(integer > 2) go to 70
      call qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24)

    // integer = 1  (or 2)
      do 30 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*ri(i)
        integeral24 = integeral24+chebyshevDegree24(i)*ri(i)
   30 continue
      do 40 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*ri(i)
   40 continue
      if(integer == 1) go to 130

    // integer = 2
      dc = alog(br-bl)
      integral = integeral24*dc
      m_estimatedError = abs((integeral24-integeral12)*dc)
      integeral12 = 0.
      integeral24 = 0.
      do 50 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rg(i)
        integeral24 = integeral12+chebyshevDegree24(i)*rg(i)
   50 continue
      do 60 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rg(i)
   60 continue
      go to 130

    // Compute the chebyshev series expansion of the following function
    // f4 = f1*log(0.5*(upperLimit+upperLimit-br-a)-0.5*(br-a)*x)
   70 fValue(1) = fValue(1)*alog(fix-halfLength)
      fValue(13) = fValue(13)*alog(fix)
      fValue(25) = fValue(25)*alog(fix+halfLength)
      do 80 i=2,12
        u = halfLength*x(i-1)
        isym = 26-i
        fValue(i) = fValue(i)*alog(fix-u)
        fValue(isym) = fValue(isym)*alog(fix+u)
   80 continue
      call qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24)

    // integer = 3  (or 4)
      do 90 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*ri(i)
        integeral24 = integeral24+chebyshevDegree24(i)*ri(i)
   90 continue
      do 100 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*ri(i)
  100 continue
      if(integer == 3) go to 130

    // integer = 4
      dc = alog(br-bl)
      integral = integeral24*dc
      m_estimatedError = abs((integeral24-integeral12)*dc)
      integeral12 = 0.
      integeral24 = 0.
      do 110 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rg(i)
        integeral24 = integeral24+chebyshevDegree24(i)*rg(i)
  110 continue
      do 120 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rg(i)
  120 continue
  130 integral = (integral+integeral24)*factor
      m_estimatedError = (m_estimatedError+abs(integeral24-integeral12))*factor
      go to 270

    // This part of the program is executed only if upperLimit = br.
    // Compute the chebyshev series expansion of the following function
    // f2 = (0.5*(upperLimit+bl-a-a)+0.5*(upperLimit-bl)*x)**alpha*f(0.5*(upperLimit-bl)*x+0.5*(upperLimit+bl))
  140 halfLength = 0.5*(br-bl)
      center = 0.5*(br+bl)
      fix = center-a
      fValue(1) = 0.5*f(halfLength+center)*(fix+halfLength)**alpha
      fValue(13) = f(center)*(fix**alpha)
      fValue(25) = 0.5*f(center-halfLength)*(fix-halfLength)**alpha
      do 150 i=2,12
        u = halfLength*x(i-1)
        isym = 26-i
        fValue(i) = f(u+center)*(fix+u)**alpha
        fValue(isym) = f(center-u)*(fix-u)**alpha
  150 continue
      factor = halfLength**(beta+1.)
      integral = 0.
      m_estimatedError = 0.
      integeral12 = 0.
      integeral24 = 0.
      if(integer == 2 || integer == 4) go to 200

    // integer = 1  (or 3)
      call qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24)
      do 160 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rj(i)
        integeral24 = integeral24+chebyshevDegree24(i)*rj(i)
  160 continue
      do 170 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rj(i)
  170 continue
      if(integer == 1) go to 260

    // integer = 3
      dc = alog(br-bl)
      integral = integeral24*dc
      m_estimatedError = abs((integeral24-integeral12)*dc)
      integeral12 = 0.
      integeral24 = 0.
      do 180 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rh(i)
        integeral24 = integeral24+chebyshevDegree24(i)*rh(i)
  180 continue
      do 190 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rh(i)
  190 continue
      go to 260

    // Compute the chebyshev series expansion of the following function
    // f3 = f2*log(0.5*(upperLimit-bl)*x+0.5*upperLimit+bl-a-a))
  200 fValue(1) = fValue(1)*alog(halfLength+fix)
      fValue(13) = fValue(13)*alog(fix)
      fValue(25) = fValue(25)*alog(fix-halfLength)
      do 210 i=2,12
        u = halfLength*x(i-1)
        isym = 26-i
        fValue(i) = fValue(i)*alog(u+fix)
        fValue(isym) = fValue(isym)*alog(fix-u)
  210 continue
      call qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24)

    // integer = 2  (or 4)
      do 220 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rj(i)
        integeral24 = integeral24+chebyshevDegree24(i)*rj(i)
  220 continue
      do 230 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rj(i)
  230 continue
      if(integer == 2) go to 260
      dc = alog(br-bl)
      integral = integeral24*dc
      m_estimatedError = abs((integeral24-integeral12)*dc)
      integeral12 = 0.
      integeral24 = 0.

    // integer = 4
      do 240 i=1,13
        integeral12 = integeral12+chebyshevDegree12(i)*rh(i)
        integeral24 = integeral24+chebyshevDegree24(i)*rh(i)
  240 continue
      do 250 i=14,25
        integeral24 = integeral24+chebyshevDegree24(i)*rh(i)
  250 continue
  260 integral = (integral+integeral24)*factor
      m_estimatedError = (m_estimatedError+abs(integeral24-integeral12))*factor
  270 return
      end



/**
* qwgts(x,lowerLimit, upperLimit,alphlowerLimit, beta,integer)
* refer to qk15w
* keywords  weight function, algebraico-logarithmic
*           end-point singularities
* 
c***purpose  this function subprogram is used together with the
c            routine qaws and defines the weight function.
*/

    Scalar alpha;
    Scalar beta;

    Scalar a;
    Scalar b;
    Scalar x;

    int integer
    
    using std::pow;
    qwgts = pow((x-a), alpha) * pow((b-x), beta);

    switch (integer)
    {

    case 1
        qwgts *= alog(xma)
        return;
    
    case 2
        qwgts *= alog(bmx)
        return;
    
    case 3
        qwgts *= alog(x-a) * alog(b-x)
        return;

    default
        return;
    }
