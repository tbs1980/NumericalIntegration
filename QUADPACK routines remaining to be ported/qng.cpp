/**
 * \file
 * \brief - The routine calculates an approximation integral to a given definite integeral i = integeral of f over (lowerLimit, upperLimit),
 *          hopefully satisfying following claim for accuracy abs(i-integral) <= max(desiredAbsoluteError,desiredRelativeError*abs(i)).
 *
 * \details - This routine is a non-adaptive, gauss-kronrod(patterson) type, automatic integeration routine for use with smooth integerands.
 *
 * \sa R. Piessens, E. de Doncker-Kapenger, C. Ueberhuber, D. Kahaner, QUADPACK, A Subroutine Package for Automatic integeration, Springer Verlag, 1983.
 *
 * \param[] The following data statements contain the abscissae and weights of the integeration rules used.
 * \param[] x1 - Abscissae common to the 10, 21, 43 and 87 point rule
 * \param[] x2 - Abscissae common to the 10, 21, 43 and 87 point rule
 * \param[] x3 - Abscissae common to the 10, 21, 43 and 87 point rule
 * \param[] x4 - Abscissae of the 87-point rule
 * \param[] w10 - Weights of the 10-point formula
 * \param[] w21a - Weights of the 21-point formula for abscissae x1
 * \param[] w21b - Weights of the 21-point formula for abscissae x2
 * \param[] w43a - Weights of the 43-point formula for abscissae x1, x3
 * \param[] w43b - Weights of the 43-point formula for abscissae x3
 * \param[] w87a - Weights of the 87-point formula for abscissae x1, x2, x3
 * \param[] w87b - Weights of the 87-point formula for abscissae x4
 * \param[] absc - abscissa
 * \param[] fValue - function value
 * \param[] savfun - array of function values which have already been computed
 * \param[] res10 - 10-point gauss integral
 * \param[] res21 - 21-point kronrod integral
 * \param[] res43 - 43-point integral
 * \param[] res87 - 87-point integral
 * \param[] absIntegral - approximation to the integeral of abs(f)
 * \param[] m_estimatedError - approximation to the integeral of abs(f-i/(upperLimit-lowerLimit))
 *
 * \returns The approximation to the integeral.
 */ 

    absc,
    center,
    dhalfLength,
    desiredAbsoluteError,
    desiredRelativeError,

    fCenter,

    fValue,
    fValue1,
    fValue2,

    halfLength,
    integral,

    res10,
    res21,
    res43,
    res87,

    absIntegral,
    m_estimatedError,

    resultKronrodh,
    r1mach,
    savfun,

    int ipx;
    int k;
    int l;

    fv1(5),
    fv2(5),
    fv3(5),
    fv4(5),

    x1(5),
    x2(5),
    x3(11),
    x4(22),

    w10(5),
    w21a(5),
    w21b(6),
    w43a(10),
    w43b(12),
    w87a(21),
    w87b(23),

    savfun(21)

    data x1
    {
        0.9739065285171717,
        0.8650633666889845,
        0.6794095682990244,
        0.4333953941292472,
        0.1488743389816312
    }.finished();

    data x2
    {
        0.9956571630258081,
        0.9301574913557082,
        0.7808177265864169,
        0.5627571346686047,
        0.2943928627014602
    }.finished();

    data x3
    {
        0.9993333609019321,
        0.9874334029080889,
        0.9548079348142663,
        0.9001486957483283,
        0.8251983149831142,
        0.7321483889893050,
        0.6228479705377252,
        0.4994795740710565,
        0.3649016613465808,
        0.2222549197766013,
        0.07465061746138332
    }.finished();

    data x4
    {
        0.9999029772627292,
        0.9979898959866787,
        0.9921754978606872,
        0.9813581635727128,
        0.9650576238583846,
        0.9431676131336706,
        0.9158064146855072,
        0.8832216577713165,
        0.8457107484624157,
        0.8035576580352310,
        0.7570057306854956,
        0.7062732097873218,
        0.6515894665011779,
        0.5932233740579611,
        0.5314936059708319,
        0.4667636230420228,
        0.3994248478592188,
        0.3298748771061883,
        0.2585035592021616,
        0.1856953965683467,
        0.1118422131799075,
        0.03735212339461987
    }.finished();

    data w10
    {
        0.06667134430868814,
        0.1494513491505806,
        0.2190863625159820,
        0.2692667193099964,
        0.2955242247147529
    }.finished();

    data w21a
    {
        0.03255816230796473,
        0.07503967481091995,
        0.1093871588022976,
        0.1347092173114733,
        0.1477391049013385
    }.finished();

    data w21b
    {
        0.01169463886737187,
        0.05475589657435200,
        0.09312545458369761,
        0.1234919762620659,
        0.1427759385770601,
        0.1494455540029169
    }.finished();

    data w43a
    {
        0.01629673428966656,
        0.03752287612086950,
        0.05469490205825544,
        0.06735541460947809,
        0.07387019963239395,
        0.005768556059769796,
        0.02737189059324884,
        0.04656082691042883,
        0.06174499520144256,
        0.07138726726869340
    }.finished();

    data w43b
    {
         0.001844477640212414,
         0.01079868958589165,
         0.02189536386779543,
         0.03259746397534569,
         0.04216313793519181,
         0.05074193960018458,
         0.05837939554261925,
         0.06474640495144589,
         0.06956619791235648,
         0.07282444147183321,
         0.07450775101417512,
         0.07472214751740301
    }.finished();
     
    data w87a
    {
        0.008148377384149173,
        0.01876143820156282,
        0.02734745105005229,
        0.03367770731163793,
        0.03693509982042791,
        0.002884872430211531,
        0.01368594602271270,
        0.02328041350288831,
        0.03087249761171336,
        0.03569363363941877,
        0.0009152833452022414,
        0.005399280219300471,
        0.01094767960111893,
        0.01629873169678734,
        0.02108156888920384,
        0.02537096976925383,
        0.02918969775647575,
        0.03237320246720279,
        0.03478309895036514,
        0.03641222073135179,
        0.03725387550304771
     }.finished();
     
    data w87b
    {
         0.000274145563762072,
         0.001807124155057943,
         0.004096869282759165,
         0.006758290051847379,
         0.900549957672201647,
         0.01232944765224485,
         0.01501044734638895,
         0.01754896798624319,
         0.01993803778644089,
         0.02219493596101229,
         0.02433914712600081,
         0.02637450541483921,
         0.02828691078877120,
         0.03005258112809270,
         0.03164675137143993,
         0.03305041341997850,
         0.03425509970422606,
         0.03526241266015668,
         0.03607698962288870,
         0.03669860449845609,
         0.03712054926983258,
         0.03733422875193504,
         0.03736107376267902
     }.finished();
     
      integral = 0.;
      m_estimatedError = 0.;
      m_numEvaluations = 0;
      m_errorCode = 6;

      if(desiredAbsoluteError <= 0. && desiredRelativeError < (std::max)(0.5e-14,5.*Eigen::NumTraits<Scalar>::epsilon()))
      {
        return;
      }
      
      halfLength = 0.5*(upperLimit-lowerLimit)
      dhalfLength = abs(halfLength)
      center = 0.5*(upperLimit+lowerLimit)
     fCenter = f(center)
      m_numEvaluations = 21
      m_errorCode = 1

    // Compute the integeral using the 10- and 21-point formula.
do 70 l = 1,3
      go to (5,25,45),l
    5 res10 = 0.
      res21 = w21b(6)* fCenter
      absIntegral = w21b(6)*abs( fCenter)
      do 10 k=1,5
        absc = halfLength*x1(k)
        fValue1 = f(center+absc)
        fValue2 = f(center-absc)
        fValue = fValue1+fValue2
        res10 = res10+w10(k)*fValue
        res21 = res21+w21a(k)*fValue
        absIntegral = absIntegral+w21a(k)*(abs(fValue1)+abs(fValue2))
        savfun(k) = fValue
        fv1(k) = fValue1
        fv2(k) = fValue2
   10 continue
      ipx = 5
      do 15 k=1,5
        ipx = ipx+1
        absc = halfLength*x2(k)
        fValue1 = f(center+absc)
        fValue2 = f(center-absc)
        fValue = fValue1+fValue2
        res21 = res21+w21b(k)*fValue
        absIntegral = absIntegral+w21b(k)*(abs(fValue1)+abs(fValue2))
        savfun(ipx) = fValue
        fv3(k) = fValue1
        fv4(k) = fValue2
   15 continue

    // Test for convergence.
      integral = res21*halfLength
      absIntegral = absIntegral*dhalfLength
      resultKronrodh = 0.5*res21
      m_estimatedError = w21b(6)*abs( fCenter-resultKronrodh)
      do 20 k = 1,5
        m_estimatedError = m_estimatedError+w21a(k)*(abs(fv1(k)-resultKronrodh)+abs(fv2(k)-resultKronrodh))+w21b(k)*(abs(fv3(k)-resultKronrodh)+abs(fv4(k)-resultKronrodh))
   20 continue
      m_estimatedError = abs((res21-res10)*halfLength)
      m_estimatedError = m_estimatedError*dhalfLength
      go to 65
    // Compute the integeral using the 43-point formula.
   25 res43 = w43b(12)* fCenter
      m_numEvaluations = 43
      do 30 k=1,10
        res43 = res43+savfun(k)*w43a(k)
   30 continue
      do 40 k=1,11
        ipx = ipx+1
        absc = halfLength*x3(k)
        fValue = f(absc+center)+f(center-absc)
        res43 = res43+fValue*w43b(k)
        savfun(ipx) = fValue
   40 continue

    // Test for convergence.
      integral = res43*halfLength
      m_estimatedError = abs((res43-res21)*halfLength)
      go to 65

    // Compute the integeral using the 87-point formula.
   45 res87 = w87b(23)* fCenter
      m_numEvaluations = 87
      do 50 k=1,21
        res87 = res87+savfun(k)*w87a(k)
   50 continue
      do 60 k=1,22
        absc = halfLength*x4(k)
        res87 = res87+w87b(k)*(f(absc+center)+f(center-absc))
   60 continue
      integral = res87*halfLength
      m_estimatedError = abs((res87-res43)*halfLength)
   65 if(m_estimatedError != 0. && m_estimatedError != 0.)
     *  m_estimatedError = m_estimatedError*(std::min)(1.,
     *  (0.2e+03*m_estimatedError/m_estimatedError)**1.5 )
      if (absIntegral > (std::numeric_limits<Scalar>::min)()/(5.*Eigen::NumTraits<Scalar>::epsilon())) m_estimatedError = (std::max)
     *  ((Eigen::NumTraits<Scalar>::epsilon()*5.)*absIntegral,m_estimatedError)
      if (m_estimatedError <= (std::max)(desiredAbsoluteError,desiredRelativeError*abs(integral))) m_errorCode = 0

    if (m_errorCode == 0)
    {
        break;
    }
    70 continue
    return;
}
