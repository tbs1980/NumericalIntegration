I compared the [GaussKronrodNonAdaptive](https://github.com/lballabio/quantlib/blob/master/QuantLib/ql/math/integrals/kronrodintegral.hpp) with our implementation as a preliminary test. According the QuantLib documentation, the method uses 15 points Gauss-Kronrod integration rule. I compred the test functions in [QuantLib test suite](https://github.com/lballabio/quantlib/blob/master/QuantLib/test-suite/integrals.cpp) and the resutls are below. The values are | expected - calculated |. The ouput is done using the following statement in C++

    std::cout<<"Normal Distribution, |calculated - expected| = "
        <<std::setprecision(10)<<std::abs(calculated - expected)<<std::endl;

| Tests         | QuantLib        | Eigen |
| ------------- |-----------------|-------|
| constant 0    | 0               | 0     |
| constant 1    | 2.997602166e-15 | 0     |
| identity      | 1.498801083e-15 | 0     |
| square        | 1.054711873e-15 | 0     |
| sine          | 5.995204333e-15 | 0     |
|Normal Distribution | 6.334248367e-05** | 2.220446049e-16 |
| Abcd2 | 6.591949209e-17 | 3.469446952e-18 |

For the Normal Distribution test the QuantLib integration limits were [-5,5] and hence the higher error. Increasing the limits beyond that gave a maximum ` maximum number of function evaluations exceeded` even when I use 10000 as the maximum number of evaluestions. The Eigen integrator does not have any issues for this function. 
