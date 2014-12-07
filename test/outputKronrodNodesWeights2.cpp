#include <NumericalIntegration.h>

#include <iostream>
#include <fstream>
#include <iomanip>

int test_values()
{
    /** 
     * \details The level of precision in the calculations requires greater precision
     *          than the number of the digits written to file to avoid roundoff errors.
     */
    int outputDigits = 256;

    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(outputDigits*2);
    typedef Eigen::QuadratureKronrod<Scalar> QuadratureKronrodValuesType;

    QuadratureKronrodValuesType::computeNodesAndWeights();

    std::ofstream fout;
    fout.open("test/testOutput/QuadratureKronrod.h");
    fout<<std::fixed<<std::setprecision(outputDigits);

    std::string gaussRule[12] = {7, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 100};
    std::string kronrodRule[12] = {15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201};

    // @TODO Create 12x4 array to hold these values 
    std::string gaussKronrodAbscissaeNames[12];
    std::string gausskronrodWeightsNames[12];
    std::string gaussAbscissaeNames[12];
    std::string gaussWeightsNames[12];

    for (size_t i=0; i<12; ++i)
    {
        gaussKronrodAbscissaeNames(i) = "abscissaeGaussKronrod";
        gaussKronrodAbscissaeNames(i) += kronrodRule(i);
        
        gaussKronrodWeightsNames(i) = "weightsGaussKronrod";
        gaussKronrodWeightsNames(i) += kronrodRule(i);
        
        gaussAbscissaeNames(i) = "abscissaeGauss";
        gaussAbscissaeNames(i) += kronrodRule(i);

        gaussWeightsNames(i) = "weightsGauss";
        gaussWeightsNames(i) += kronrodRule(i);
    }

//----------------------------------Begin File Header Information--------------------------------//

    fout << "/**\n * \\file QuadratureKronrod.h\n"
         << " * \\sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package\n"
         << " *     for Automatic Integration, Springer Verlag, 1983.\n"
         << " */\n\n#ifndef EIGEN_QUADRATURE_KRONROD_H\n#define EIGEN_QUADRATURE_KRONROD_H\n\n"
         << "namespace Eigen\n{\n\n/**\n * \\brief The abscissae and weights are given for the interval (-1,1).\n"
         << " *        Because of symmetry, only the positive abscissae and their\n"
         << " *        corresponding weights are given.\n*\n";

    for (size_t i=0; i<12; ++1)
    {
        // " * \\param abscissaeGaussKronrodXX  The abscissae of the X point kronrod rule."
        fout << " * \\param abscissaeGaussKronrod" << kronrodRule(i) << "  The abscissae of the " << kronrodRule(i) << " point kronrod rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++1)
    {
        // " * \\param weightsGaussKronrodXX  The weights of the X point kronrod rule."
        fout << " * \\param weightsGaussKronrod" << kronrodRule(i) << "  The weights of the " << kronrodRule(i) << " point kronrod rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++1)
    {
        // " * \\param abscissaeGaussXX  The abscissae of the X point gauss rule."
        fout << " * \\param abscissaeGauss" << kronrodRule(i) << "  The abscissae of the " << gaussRule(i) << " point gauss rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++1)
    {
        // " * \\param weightsGaussXX  The weights of the X point gauss rule."
        fout << " * \\param weightsGauss" << kronrodRule(i) << "  The weights of the " << gaussRule(i) << " point gauss rule.\n";
    }

    fout << " *\ntemplate <typename Scalar>\nclass QuadratureKronrod\n{\npublic:\n";
    
    for (size_t i=0; i<12; ++1)
    {
        fout << "\tstatic Array<Scalar, " << kronrodRule(i) << ", 1> " << gaussKronrodAbscissaeNames(i) << std::endl;
             << "\tstatic Array<Scalar, " << kronrodRule(i) << ", 1> " << gausskronrodWeightsNames(i) << std::endl;
             << "\tstatic Array<Scalar, " << gaussRule(i) << ", 1> " << gaussAbscissaeNames(i) << std::endl;
             << "\tstatic Array<Scalar, " << gaussRule(i) << ", 1> " << gaussWeightsNames(i) << std::endl;
    }

    fout << "\ttypedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;\n"
         << "\ttypedef Kronrod::Piessens<Scalar> PiessensPolicy;\n"
         << "\ttypedef typename LaurieGautschiPolicy::VectorType VectorType;\n\n"
         << "\tstatic bool compute;\n\n";

    for (size_t i=0; i<12; ++1)
    {
        fout << "\tstatic void computeNodesAndWeights()\n\t\t{\n\t\tif(compute)\n\t\t{\n"
             << "\t\t\tQuadratureKronrod::computeForRule<" << gaussRule(i) << ">(" 
             << gaussKronrodAbscissaeNames(i) << ", " << gausskronrodWeightsNames(i) << ", "
             << gaussAbscissaeNames(i) << ", " << gaussWeightsNames(i) << ");\n";
    }

    fout << "\n\t\t\tcompute = false;\n\t\t}\n\t}\n\n"
         << "\ttemplate <int N>\n"
         << "\tstatic void computeForRule(Array<Scalar, N+1, 1>& kronrodAbscissae, Array<Scalar, N+1, 1>& kronrodWeights,\n"
         << "\t\t\t\t\t\t\t   Array<Scalar, (N+1)/2, 1>& gaussAbscissae, Array<Scalar, (N+1)/2, 1>& gaussWeights)\n\t\t{\n"
         << "\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> xGK;\n\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> wGK;\n"
         << "\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> xG;\n\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> wG;\n\n"
         << "\t\tLaurieGautschiPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n"
         << "\t\t//PiessensPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n\n"
         << "\t\tfor(size_t i=0; i<N+1; ++i)\n\t\t{\n\t\t\tkronrodAbscissae(i) = xGK(i);\n\t\t\tkronrodWeights(i) =  wGK(i);\n\t\t}\n\n"
         << "\t\tfor(size_t i=0; i<(N+1)/2; ++i)\n\t\t{\n\t\t\tgaussAbscissae(i) = xG(i);\n\t\t\tgaussWeights(i) = wG(i);\n\t\t}\n\t}\n};\n\n"
         << "template <typename Scalar>\nbool QuadratureKronrod<Scalar>::compute = true;\n\n";

//-----------------------------------End File Header Information---------------------------------//

    QuadratureKronrodValuesType::abscissaeGaussKronrod15(i)

    Eigen::Array<Scalar,Eigen::Dynamic,1>  gaussKronrodAbscissae;
    Eigen::Array<Scalar,Eigen::Dynamic,1> kronrodWeights;
    Eigen::Array<Scalar,Eigen::Dynamic,1>  gaussAbscissae;
    Eigen::Array<Scalar,Eigen::Dynamic,1> gaussWeights;

    int N;

    for (size_t i=0; i<12; ++i)
    {
        N = gaussRule(i);
        int kronrodRule = 2 * gaussRule(i) + 1;

        for (int j=0; j<2; ++j)
        {
            if (j==0)
            {
                fout<< "template <typename Scalar>\n"
                    << "Array<Scalar," << N+1 << ", 1> QuadratureKronrod<Scalar>::" << gaussKronrodAbscissaeNames(i) << kronrodRule << " =\n"
                    << "  (Array<Scalar," << N+1 << ", 1>() <<\n";
            }
            if (j==1)
            {
                fout<<"template <typename Scalar>\n"
                    << "Array<Scalar," << N+1 << ", 1> QuadratureKronrod<Scalar>::" << gaussKronrodWeightsNames(i) << kronrodRule << " =\n"
                    << "  (Array<Scalar," << N+1 << ", 1>() <<\n";
            }

            for(int i=0;i<N+1;++i)
            {
                // @TODO The next line needs to be replaced with functional code
                fout<<"\t"<<QuadratureKronrodValuesType::abscissaeGaussKronrod15(i);
                
                if(i==N)
                {
                    fout<<"\n  ).finished();\n\n";
                    break;
                }

                fout<<",\n";
            }
        }

        for (int j=2; j<4; ++j)
        {
            if (j==2)
            {
                fout<<"template <typename Scalar>\n"
                    << "Array<Scalar," << (N+1)/2 << ", 1> QuadratureKronrod<Scalar>::"<< gaussAbscissaeNames(i) << kronrodRule << " =\n"
                    << "  (Array<Scalar," << (N+1)/2 << ", 1>() <<\n";
            }
            if (j==3)
            {
                fout<<"template <typename Scalar>\n"
                    << "Array<Scalar," << (N+1)/2 << ", 1> QuadratureKronrod<Scalar>::" << gaussWeightsNames(i) << kronrodRule << " =\n"
                    << "  (Array<Scalar," << (N+1)/2 << ", 1>() <<\n";
            }

            for(int i=0;i<(N+1)/2;++i)
            {
                // @TODO The next line needs to be replaced with functional code
                fout<<"\t"<<QuadratureKronrodValuesType::abscissaeGauss15(i);
                
                if(i==(N+1)/2 - 1)
                {
                    fout<<"\n  ).finished();\n\n";
                    break;
                }

                fout<<",\n";
            }
        }
    }

    fout.close();
    std::cout << "\n  Kronrod Nodes and Weights written to file \"test/testOutput/QuadratureKronrod.h\"\n\n.";

    return EXIT_SUCCESS;
}

int main(void)
{
    int ret=EXIT_SUCCESS;
    test_values();
    return ret;
}
