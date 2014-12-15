/**
* subroutine qcheb(x,fValue,chebyshevDegree12,chebyshevDegree24)
*
* keywords - Chebyshev series expansion, fast fourm_errorCode transform
*
* purpose - This routine computes the chebyshev series expansion of degrees 12 and 24 of a function using a
*            fast fourm_errorCode transform method f(x) = sum(k=1,..,13) (chebyshevDegree12(k)*t(k-1,x)),
*            f(x) = sum(k=1,..,25) (chebyshevDegree24(k)*t(k-1,x)), where t(k,x) is the chebyshev polynomial of degree k.
*
* description - Chebyshev series expansion
*
*   x      - Vector of dimension 11 containing the values cos(k*M_PI/24), from k = 1 to 11
*
*   fValue - Vector of dimension 25 containing the function values at the points
*            (upperLimit+lowerLimit+(upperLimit-lowerLimit)*cos(k*M_PI/24))/2, from k = 0 to 24,
*            where (lowerLimit, upperLimit) is the approximation interval. fValue(1) and fValue(25)
*            are divided by two (these values are destroyed at output).
*
*   chebyshevDegree12 - Vector of dimension 13 containing the Chebyshev coefficients for degree 12
*
*   chebyshevDegree24 - Vector of dimension 25 containing the Chebyshev coefficients for degree 24
*/
    chebyshevDegree12(13);
    chebyshevDegree24(25);

    fValue(25);
    v(12);
    x(11);
    
    Scalar alam;
    Scalar alam1;
    Scalar alam2;

    Scalar part1;
    Scalar part2;
    Scalar part3;
    
    int i
    int j;

    for (size_t i=0; i<12; ++i)
    {
        j = 26-i;
        v(i) = fValue(i) - fValue(j);
        fValue(i) = fValue(i) + fValue(j);
    }

    alam1 = v(1) - v(9);
    alam2 = x(6)*(v(3) - v(7) - v(11));

    chebyshevDegree12(4) = alam1 + alam2;
    chebyshevDegree12(10) = alam1 - alam2;

    alam1 = v(2) - v(8) - v(10);
    alam2 = v(4) - v(6) - v(12);
    alam = x(3)*alam1 + x(9)*alam2;

    chebyshevDegree24(4) = chebyshevDegree12(4) + alam;
    chebyshevDegree24(22) = chebyshevDegree12(4) - alam;

    alam = x(9)*alam1 - x(3)*alam2;

    chebyshevDegree24(10) = chebyshevDegree12(10)+alam;
    chebyshevDegree24(16) = chebyshevDegree12(10)-alam;

    part1 = x(4)*v(5);
    part2 = x(8)*v(9);
    part3 = x(6)*v(7);

    alam1 = v(1) + part1 + part2;
    alam2 = x(2)*v(3) + part3 + x(10)*v(11);

    chebyshevDegree12(2) = alam1 + alam2;
    chebyshevDegree12(12) = alam1 - alam2;
    
    alam = x(1)*v(2) + x(3)*v(4) + x(5)*v(6) + x(7)*v(8) + x(9)*v(10) + x(11)*v(12);
    
    chebyshevDegree24(2) = chebyshevDegree12(2) + alam;
    chebyshevDegree24(24) = chebyshevDegree12(2) - alam;

    alam = x(11)*v(2) - x(9)*v(4) + x(7)*v(6) - x(5)*v(8) + x(3)*v(10) - x(1)*v(12);
    
    chebyshevDegree24(12) = chebyshevDegree12(12) + alam;
    chebyshevDegree24(14) = chebyshevDegree12(12) - alam;
    
    alam1 = v(1) - part1+part2;
    alam2 = x(10)*v(3) - part3 + x(2)*v(11);
    
    chebyshevDegree12(6) = alam1 + alam2;
    chebyshevDegree12(8) = alam1 - alam2;
    
    alam = x(5)*v(2) - x(9)*v(4) - x(1)*v(6) - x(11)*v(8) + x(3)*v(10) + x(7)*v(12);

    chebyshevDegree24(6) = chebyshevDegree12(6) + alam;
    chebyshevDegree24(20) = chebyshevDegree12(6) - alam;
    
    alam = x(7)*v(2) - x(3)*v(4) - x(11)*v(6) + x(1)*v(8) - x(9)*v(10) - x(5)*v(12);

    chebyshevDegree24(8) = chebyshevDegree12(8) + alam;
    chebyshevDegree24(18) = chebyshevDegree12(8) - alam;
    
    for(size_t i=0; i<6; ++i)
    {
        j = 14 - i;
        v(i) = fValue(i) - fValue(j);
        fValue(i) = fValue(i) + fValue(j);
    }

    alam1 = v(1) + x(8)*v(5);
    alam2 = x(4)*v(3);

    chebyshevDegree12(3) = alam1 + alam2;
    chebyshevDegree12(7) = v(1) - v(5);
    chebyshevDegree12(11) = alam1 - alam2;


    alam = x(2)*v(2) + x(6)*v(4) + x(10)*v(6);

    chebyshevDegree24(3) = chebyshevDegree12(3) + alam;
    chebyshevDegree24(23) = chebyshevDegree12(3) - alam;

    alam = x(6)*(v(2) - v(4) - v(6));

    chebyshevDegree24(7) = chebyshevDegree12(7) + alam;
    chebyshevDegree24(19) = chebyshevDegree12(7) - alam;

    alam = x(10)*v(2) - x(6)*v(4) + x(2)*v(6);

    chebyshevDegree24(11) = chebyshevDegree12(11) + alam;
    chebyshevDegree24(15) = chebyshevDegree12(11) - alam;
    
    for (size_t i=0; i<3; ++i)
    {
        j = 8-i
        v(i) = fValue(i) - fValue(j);
        fValue(i) = fValue(i) + fValue(j);
    }

    chebyshevDegree12(5) = v(1) + x(8)*v(3);
    chebyshevDegree12(9) = fValue(1) - x(8)*fValue(3);

    alam = x(4)*v(2);

    chebyshevDegree24(5) = chebyshevDegree12(5) + alam;
    chebyshevDegree24(21) = chebyshevDegree12(5) - alam;

    alam = x(8)*fValue(2) - fValue(4);

    chebyshevDegree12(1) = fValue(1) + fValue(3);
    chebyshevDegree12(13) = v(1) - v(3);

    chebyshevDegree24(9) = chebyshevDegree12(9) + alam;
    chebyshevDegree24(17) = chebyshevDegree12(9) - alam;

    alam = fValue(2) + fValue(4);

    chebyshevDegree24(1) = chebyshevDegree12(1) + alam;
    chebyshevDegree24(13) = chebyshevDegree12(13);
    chebyshevDegree24(25) = chebyshevDegree12(1) - alam;

    alam = 1./6.;
    
    for (size_t i=1; i<12; ++i)
    {
        chebyshevDegree12(i) = chebyshevDegree12(i)*alam;
    }

    alam = 0.5*alam;
    chebyshevDegree12(1) = chebyshevDegree12(1)*alam;
    chebyshevDegree12(13) = chebyshevDegree12(13)*alam;
    
    for (size_t i=1; i<24; ++i)
    {
        chebyshevDegree24(i) = chebyshevDegree24(i)*alam;
    }

    chebyshevDegree24(1) = 0.5 * alam * chebyshevDegree24(1);
    chebyshevDegree24(25) = 0.5 * alam * chebyshevDegree24(25);
    return;
}
