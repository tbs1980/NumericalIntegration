#include <NIHeaders.h>

#include <iostream>
#include <fstream>
#include <iomanip>

int test_values()
{
    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(256);
    typedef Eigen::QuadratureKronrod<Scalar> QuadratureKronrodValuesType;

    QuadratureKronrodValuesType::ComputeNodesAndWeights();

    int outputDigits = 256;
    ofstream fout;
    fout.open("KronrodNodesAndWeights.txt");

    //--------15--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 8, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod15 ="
        <<"\n  (Array<Scalar, 8, 1>() <<\n";

    for(size_t i=0;i<8;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod15(i);
        
        if(i<7)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 8, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod15 ="
        <<"\n  (Array<Scalar, 8, 1>() <<\n";

    for(size_t i=0;i<8;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod15(i);
        
        if(i<7)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 4, 1> QuadratureKronrod<Scalar>::weightsGauss15 ="
        <<"\n  (Array<Scalar, 4, 1>() <<\n";

    for(size_t i=0;i<4;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss15(i);
        
        if(i<3)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------21--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 11, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod21 ="
        <<"\n  (Array<Scalar, 11, 1>() <<\n";

    for(size_t i=0;i<11;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod21(i);

        if(i<10)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }
    
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 11, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod21 ="
        <<"\n  (Array<Scalar, 11, 1>() <<\n";

    for(size_t i=0;i<11;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod21(i);

        if(i<10)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 5, 1> QuadratureKronrod<Scalar>::weightsGauss21 ="
        <<"\n  (Array<Scalar, 5, 1>() <<\n";

    for(size_t i=0;i<5;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss21(i);

        if(i<4)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------31--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 16, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod31 ="
        <<"\n  (Array<Scalar, 16, 1>() <<\n";

    for(size_t i=0;i<16;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod31(i);

        if(i<15)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 16, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod31 ="
        <<"\n  (Array<Scalar, 16, 1>() <<\n";

    for(size_t i=0;i<16;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod31(i);

        if(i<15)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 8, 1> QuadratureKronrod<Scalar>::weightsGauss31 ="
        <<"\n  (Array<Scalar, 8, 1>() <<\n";

    for(size_t i=0;i<8;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss31(i);

        if(i<7)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------41--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 21, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod41 ="
        <<"\n  (Array<Scalar, 21, 1>() <<\n";

    for(size_t i=0;i<21;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod41(i);

        if(i<20)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 21, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod41 ="
        <<"\n  (Array<Scalar, 21, 1>() <<\n";

    for(size_t i=0;i<21;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod41(i);

        if(i<20)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 10, 1> QuadratureKronrod<Scalar>::weightsGauss41 ="
        <<"\n  (Array<Scalar, 10, 1>()  <<\n";

    for(size_t i=0;i<10;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss41(i);

        if(i<9)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------51--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 26, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod51 ="
        <<"\n  (Array<Scalar, 26, 1>() <<\n";
    for(size_t i=0;i<26;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod51(i);

        if(i<25)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 26, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod51 ="
        <<"\n  (Array<Scalar, 26, 1>() <<\n";

    for(size_t i=0;i<26;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod51(i);

        if(i<25)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 13, 1> QuadratureKronrod<Scalar>::weightsGauss51 ="
        <<"\n  (Array<Scalar, 13, 1>() <<\n";

    for(size_t i=0;i<13;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss51(i);

        if(i<12)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------61--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 31, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod61 ="
        <<"\n  (Array<Scalar, 31, 1>() <<\n";

    for(size_t i=0;i<31;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod61(i);

        if(i<30)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 31, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod61 ="
        <<"\n  (Array<Scalar, 31, 1>() <<\n";
    
    for(size_t i=0;i<31;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod61(i);

        if(i<30)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 15, 1> QuadratureKronrod<Scalar>::weightsGauss61 ="
        <<"\n  (Array<Scalar, 15, 1>() <<\n";

    for(size_t i=0;i<15;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss61(i);

        if(i<14)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------71--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 36, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod71 ="
        <<"\n  (Array<Scalar, 36, 1>() <<\n";

    for(size_t i=0;i<36;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod71(i);

        if(i<35)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 36, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod71 ="
        <<"\n  (Array<Scalar, 36, 1>() <<\n";
    
    for(size_t i=0;i<36;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod71(i);

        if(i<35)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 18, 1> QuadratureKronrod<Scalar>::weightsGauss71 ="
        <<"\n  (Array<Scalar, 18, 1>() <<\n";

    for(size_t i=0;i<18;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss71(i);

        if(i<17)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------81--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 41, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod81 ="
        <<"\n  (Array<Scalar, 41, 1>() <<\n";
    for(size_t i=0;i<41;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod81(i);

        if(i<40)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 41, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod81 ="
        <<"\n  (Array<Scalar, 41, 1>() <<\n";
    
    for(size_t i=0;i<41;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod81(i);

        if(i<40)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 20, 1> QuadratureKronrod<Scalar>::weightsGauss81 ="
        <<"\n  (Array<Scalar, 20, 1>() <<\n";

    for(size_t i=0;i<20;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss81(i);

        if(i<19)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }
    
    //--------91--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 46, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod91 ="
        <<"\n  (Array<Scalar, 46, 1>() <<\n";

    for(size_t i=0;i<46;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod91(i);

        if(i<45)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 46, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod91 ="
        <<"\n  (Array<Scalar, 46, 1>() <<\n";
    
    for(size_t i=0;i<46;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod91(i);

        if(i<45)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 23, 1> QuadratureKronrod<Scalar>::weightsGauss91 ="
        <<"\n  (Array<Scalar, 23, 1>() <<\n";

    for(size_t i=0;i<23;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss91(i);

        if(i<22)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------101--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 51, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod101 ="
        <<"\n  (Array<Scalar, 51, 1>() <<\n";

    for(size_t i=0;i<51;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod101(i);

        if(i<50)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 51, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod101 ="
        <<"\n  (Array<Scalar, 51, 1>() <<\n";
    
    for(size_t i=0;i<51;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod101(i);

        if(i<50)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 25, 1> QuadratureKronrod<Scalar>::weightsGauss101 ="
        <<"\n  (Array<Scalar, 25, 1>() <<\n";

    for(size_t i=0;i<25;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss101(i);

        if(i<24)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------121--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 61, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod121 ="
        <<"\n  (Array<Scalar, 61, 1>() <<\n";

    for(size_t i=0;i<61;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod121(i);

        if(i<60)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 61, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod121 ="
        <<"\n  (Array<Scalar, 61, 1>() <<\n";
    
    for(size_t i=0;i<61;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod121(i);

        if(i<60)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }
    
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 30, 1> QuadratureKronrod<Scalar>::weightsGauss121 ="
        <<"\n  (Array<Scalar, 30, 1>() <<\n";

    for(size_t i=0;i<30;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss121(i);

        if(i<29)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    //--------201--------//
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 101, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod201 ="
        <<"\n  (Array<Scalar, 101, 1>() <<\n";

    for(size_t i=0;i<101;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod201(i);

        if(i<100)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 101, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod201 ="
        <<"\n  (Array<Scalar, 101, 1>() <<\n";
    
    for(size_t i=0;i<101;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGaussKronrod201(i);

        if(i<100)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }
    
    fout<<"template <typename Scalar>"
        <<"\nArray<Scalar, 50, 1> QuadratureKronrod<Scalar>::weightsGauss201 ="
        <<"\n  (Array<Scalar, 50, 1>() <<\n";

    for(size_t i=0;i<50;++i)
    {
        fout<<std::setprecision(outputDigits)<<"\t"<<QuadratureKronrodValuesType::weightsGauss201(i);

        if(i<49)
        {
            fout<<",\n";
        }
        else
        {
            fout<<"\n  ).finished();\n\n";
        }
    }

    fout.close();
    std::cout<<std::endl<<"  Kronrod Nodes and Weights written to file \"KronrodNodesAndWeights.txt\"\n\n.";

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=EXIT_SUCCESS;
    test_values();
    return ret;
}
