#include <NIHeaders.hpp>
#include <iostream>
#include <iomanip>

int compare_codes(void)
{
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    //RealType::set_default_prec(128);
    RealType::set_default_prec(256);

    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;
    typedef LaurieGautschiPolicy::VectorType VectorType;

    const IndexType N=7;
    const int outputIntegers = 33;

    Eigen::Array<RealType,Dynamic,2> ans;
    ans = Kronrod::multiPrecisionKronrod<RealType>(N);

    std::cout<<std::fixed;
    std::cout<<std::endl<<"MS Laurie Gautschi"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << ans(i,0)
		  << "\t" << ans(i,1) << std::endl;
    }

    VectorType x=VectorType::Zero(2*N+1);
    VectorType w=VectorType::Zero(2*N+1);
    LaurieGautschiPolicy::mpkonrad(N,x,w);

    std::cout<<"\nSTB Laurie Gautschi"<<std::endl;
    for(IndexType i=0;i<x.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << x(i)
		  << "\t" << w(i) << std::endl;
    }

    Eigen::Array<RealType, Dynamic, 1> xGKPosAndNeg
	= Eigen::Array<RealType, Dynamic, 1>::Zero(x.rows());
    Eigen::Array<RealType, Dynamic, 1> wGKPosAndNeg
	= Eigen::Array<RealType, Dynamic, 1>::Zero(x.rows());

    Eigen::Array<RealType, Dynamic, 1> xGK;
    Eigen::Array<RealType, Dynamic, 1> wGK;
    Eigen::Array<RealType, Dynamic, 1> wG;

    Kronrod::kronrod(N, xGK,  wGK, wG);

    for(IndexType i=0; i<xGK.rows(); ++i)
    {
        xGKPosAndNeg(i) = -xGK(i);
        wGKPosAndNeg(i) = wGK(i);
    }

    for(IndexType i=0; i<xGK.rows(); ++i)
    {
        xGKPosAndNeg(xGKPosAndNeg.rows()-1-i) = xGK(i);
        wGKPosAndNeg(wGKPosAndNeg.rows()-1-i) = wGK(i);
    }

    std::cout << "\nMS Piessens" << std::endl;
    for(IndexType i=0;i<xGKPosAndNeg.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << xGKPosAndNeg(i)
		  << "\t"<<wGKPosAndNeg(i) << std::endl;
    }

    std::cout << std::endl;

    for(IndexType i=0;i<wG.rows();++i)
    {
        std::cout << std::setprecision(outputIntegers) << wG(i) << std::endl;
    }



    std::cout<<std::endl;

/*------------------------Output of Differences Between Approaches-------------------------------//

    std::cout<<"\nSolution Differences: MS - STB Laurie Gautschi"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout << std::setprecision(15) << ans(i,0) - x(i)
		  << "\t" << ans(i,1) - w(i) << std::endl;
    }

    std::cout<<"\nSolution Differences: MS Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<ans.rows();++i)
    {
        std::cout << std::setprecision(15) << ans(i,0) - xGKPosAndNeg(i)
		  << "\t" << ans(i,1) - wGKPosAndNeg(i) << std::endl;
    }

    std::cout<<"\nSolution Differences: STB Laurie Gautschi - MS Piessens"<<std::endl;
    for(int i=0;i<x.rows();++i)
    {
        std::cout  << std::setprecision(15) << x(i) - xGKPosAndNeg(i)
		   << "\t" << w(i) - wGKPosAndNeg(i) << std::endl;
    }

//-----------------------End Output of Differences Between Approaches----------------------------*/

    return EXIT_SUCCESS;
}

int quad_interface_15(void)
{
    std::cout<<"\n GaussKronrod15\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,8,1> abscissaeGaussKronrod15;
    Eigen::Array<RealType,8,1> weightsGaussKronrod15;

    LaurieGautschiPolicy::mpkonrad15(abscissaeGaussKronrod15,weightsGaussKronrod15);

    std::cout<<std::fixed;
    for(IndexType i=0;i<8;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod15(i)<<"\t"<<weightsGaussKronrod15(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}

int quad_interface_21(void)
{
    std::cout<<"\n GaussKronrod21\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,11,1> abscissaeGaussKronrod21;
    Eigen::Array<RealType,11,1> weightsGaussKronrod21;

    LaurieGautschiPolicy::mpkonrad21(abscissaeGaussKronrod21,weightsGaussKronrod21);

    std::cout<<std::fixed;
    for(IndexType i=0;i<11;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod21(i)<<"\t"<<weightsGaussKronrod21(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}

int quad_interface_31(void)
{
    std::cout<<"\n GaussKronrod31\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,16,1> abscissaeGaussKronrod31;
    Eigen::Array<RealType,16,1> weightsGaussKronrod31;

    LaurieGautschiPolicy::mpkonrad31(abscissaeGaussKronrod31,weightsGaussKronrod31);

    std::cout<<std::fixed;
    for(IndexType i=0;i<16;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod31(i)<<"\t"<<weightsGaussKronrod31(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}

int quad_interface_41(void)
{
    std::cout<<"\n GaussKronrod41\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,21,1> abscissaeGaussKronrod41;
    Eigen::Array<RealType,21,1> weightsGaussKronrod41;

    LaurieGautschiPolicy::mpkonrad41(abscissaeGaussKronrod41,weightsGaussKronrod41);

    std::cout<<std::fixed;
    for(IndexType i=0;i<21;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod41(i)<<"\t"<<weightsGaussKronrod41(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}

int quad_interface_51(void)
{
    std::cout<<"\n GaussKronrod51\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,26,1> abscissaeGaussKronrod51;
    Eigen::Array<RealType,26,1> weightsGaussKronrod51;

    LaurieGautschiPolicy::mpkonrad51(abscissaeGaussKronrod51,weightsGaussKronrod51);

    std::cout<<std::fixed;
    for(IndexType i=0;i<26;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod51(i)<<"\t"<<weightsGaussKronrod51(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}

int quad_interface_61(void)
{
    std::cout<<"\n GaussKronrod61\n"<<std::endl;
    //typedef double RealType;
    typedef mpfr::mpreal RealType;
    RealType::set_default_prec(256);
    typedef Kronrod::LaurieGautschi<RealType> LaurieGautschiPolicy;
    typedef LaurieGautschiPolicy::IndexType IndexType;

    Eigen::Array<RealType,31,1> abscissaeGaussKronrod61;
    Eigen::Array<RealType,31,1> weightsGaussKronrod61;

    LaurieGautschiPolicy::mpkonrad61(abscissaeGaussKronrod61,weightsGaussKronrod61);

    std::cout<<std::fixed;
    for(IndexType i=0;i<31;++i)
    {
        std::cout << std::setprecision(33) <<abscissaeGaussKronrod61(i)<<"\t"<<weightsGaussKronrod61(i)<<std::endl;
    }

    return EXIT_SUCCESS;

}


int main(void)
{
    int ret = 0;
    //ret += compare_codes();
    ret += quad_interface_15();
    ret += quad_interface_21();
    ret += quad_interface_31();
    ret += quad_interface_41();
    ret += quad_interface_51();
    ret += quad_interface_61();

    return ret;
}
