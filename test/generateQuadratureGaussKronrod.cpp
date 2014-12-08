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

    // Set the number of output digits
    int outputDigits = 80;
    // Set this flag to 1 for LaurieGautschi Polity, or 0 for PiessensPolicy;
    bool laurieGautschiOrPiessensPolicyFlag = 1;
    
    //typedef float Scalar;
    //typedef double Scalar;
    //typedef long double Scalar;
    typedef mpfr::mpreal Scalar;
    Scalar::set_default_prec(outputDigits*2);
    typedef Eigen::QuadratureKronrod<Scalar> QuadratureKronrodValuesType;

    QuadratureKronrodValuesType::computeNodesAndWeights();

    std::ofstream fout;
    fout.open("test/testOutput/QuadratureGaussKronrod.h");
    fout<<std::fixed<<std::setprecision(outputDigits);

    int gaussRule[12] = {7, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 100};
    int kronrodRule[12] = {15, 21, 31, 41, 51, 61, 71, 81, 91, 101, 121, 201};

    // @TODO Create 12x4 array to hold these values 
    std::string gaussKronrodAbscissaeNames[12];
    std::string gaussKronrodWeightsNames[12];
    std::string gaussAbscissaeNames[12];
    std::string gaussWeightsNames[12];

    for (size_t i=0; i<12; ++i)
    {
        char ruleStr[5];
        sprintf(ruleStr, "%d", kronrodRule[i]);
        gaussKronrodAbscissaeNames[i] = "abscissaeGaussKronrod";
        gaussKronrodAbscissaeNames[i] += ruleStr;
        
        gaussKronrodWeightsNames[i] = "weightsGaussKronrod";
        gaussKronrodWeightsNames[i] +=ruleStr;
        
        gaussAbscissaeNames[i] = "abscissaeGauss";
        gaussAbscissaeNames[i] += ruleStr;

        gaussWeightsNames[i] = "weightsGauss";
        gaussWeightsNames[i] += ruleStr;
    }

//----------------------------------Begin File Header Information--------------------------------//

    fout << "/**\n * \\file QuadratureKronrod.h\n"
         << " * \\sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package\n"
         << " *     for Automatic Integration, Springer Verlag, 1983.\n"
         << " */\n\n#ifndef EIGEN_QUADRATURE_KRONROD_H\n#define EIGEN_QUADRATURE_KRONROD_H\n\n"
         << "namespace Eigen\n{\n\n/**\n * \\brief The abscissae and weights are given for the interval (-1,1).\n"
         << " *        Because of symmetry, only the positive abscissae and their\n"
         << " *        corresponding weights are given.\n *\n";

    for (size_t i=0; i<12; ++i)
    {
        // " * \\param abscissaeGaussKronrodXX  The abscissae of the X point kronrod rule."
        fout << " * \\param abscissaeGaussKronrod" << kronrodRule[i] << "  The abscissae of the " << kronrodRule[i] << " point kronrod rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++i)
    {
        // " * \\param weightsGaussKronrodXX  The weights of the X point kronrod rule."
        fout << " * \\param weightsGaussKronrod" << kronrodRule[i] << "  The weights of the " << kronrodRule[i] << " point kronrod rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++i)
    {
        // " * \\param abscissaeGaussXX  The abscissae of the X point gauss rule."
        fout << " * \\param abscissaeGauss" << kronrodRule[i] << "  The abscissae of the " << gaussRule[i] << " point gauss rule.\n";
    }
    fout << " *\n";

    for (size_t i=0; i<12; ++i)
    {
        // " * \\param weightsGaussXX  The weights of the X point gauss rule."
        fout << " * \\param weightsGauss" << kronrodRule[i] << "  The weights of the " << gaussRule[i] << " point gauss rule.\n";
    }

    fout << " */\n\ntemplate <typename Scalar>\nclass QuadratureKronrod\n{\npublic:\n";
    
    for (size_t i=0; i<12; ++i)
    {
        fout << "\tstatic Array<Scalar, " << kronrodRule[i]/2 + 1 << ", 1> " << gaussKronrodAbscissaeNames[i] << ";\n"
             << "\tstatic Array<Scalar, " << kronrodRule[i]/2 + 1 << ", 1> " << gaussKronrodWeightsNames[i] << ";\n"
             << "\tstatic Array<Scalar, " << (gaussRule[i]+1)/2 << ", 1> " << gaussAbscissaeNames[i] << ";\n"
             << "\tstatic Array<Scalar, " << (gaussRule[i]+1)/2 << ", 1> " << gaussWeightsNames[i] << ";\n\n";
    }

    fout << "\ttypedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;\n"
         << "\ttypedef Kronrod::Piessens<Scalar> PiessensPolicy;\n"
         << "\ttypedef typename LaurieGautschiPolicy::VectorType VectorType;\n\n"
         << "\tstatic bool compute;\n\n";


    fout << "\tstatic void computeNodesAndWeights()\n\t{\n\t\tif(compute)\n\t\t{\n";
    for (size_t i=0; i<12; ++i)
    {
        fout << "\t\t\tQuadratureKronrod::computeForRule<" << gaussRule[i] << ">(" 
             << gaussKronrodAbscissaeNames[i] << ", " << gaussKronrodWeightsNames[i] << ", "
             << gaussAbscissaeNames[i] << ", " << gaussWeightsNames[i] << ");\n";
    }

    fout << "\n\t\t\tcompute = false;\n\t\t}\n\t}\n\n"
         << "\ttemplate <int N>\n"
         << "\tstatic void computeForRule(Array<Scalar, N+1, 1>& kronrodAbscissae, Array<Scalar, N+1, 1>& kronrodWeights,\n"
         << "\t\t\t\t\t\t\t   Array<Scalar, (N+1)/2, 1>& gaussAbscissae, Array<Scalar, (N+1)/2, 1>& gaussWeights)\n\t{\n"
         << "\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> xGK;\n\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> wGK;\n"
         << "\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> xG;\n\t\tEigen::Array<Scalar, Eigen::Dynamic, 1> wG;\n\n";
    
    if(laurieGautschiOrPiessensPolicyFlag)
    {
        fout << "\t\tLaurieGautschiPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n"
             << "\t\t//PiessensPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n\n";
    }else
    {
        fout << "\t\t//LaurieGautschiPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n"
             << "\t\tPiessensPolicy::computeAbscissaeAndWeights((unsigned int)N,xGK,wGK,xG,wG);\n\n";
    }    
    
    fout << "\t\tfor(size_t i=0; i<N+1; ++i)\n\t\t{\n\t\t\tkronrodAbscissae(i) = xGK(i);\n\t\t\tkronrodWeights(i) =  wGK(i);\n\t\t}\n\n"
         << "\t\tfor(size_t i=0; i<(N+1)/2; ++i)\n\t\t{\n\t\t\tgaussAbscissae(i) = xG(i);\n\t\t\tgaussWeights(i) = wG(i);\n\t\t}\n\t}\n};\n\n"
         << "template <typename Scalar>\nbool QuadratureKronrod<Scalar>::compute = true;\n\n";

//-----------------------------------End File Header Information---------------------------------//


    for (size_t i=0; i<12; ++i)
    {
        int N = gaussRule[i];
        int kronrodRule = 2 * N + 1;
        Eigen::Array<Scalar,Eigen::Dynamic,1> gaussKronrodAbscissae;
        Eigen::Array<Scalar,Eigen::Dynamic,1> kronrodWeights;
        Eigen::Array<Scalar,Eigen::Dynamic,1> gaussAbscissae;
        Eigen::Array<Scalar,Eigen::Dynamic,1> gaussWeights;
        
        switch(kronrodRule)
        {
            case 15:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod15;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod15;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss15;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss15;
                break;
            case 21:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod21;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod21;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss21;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss21;
                break;
            case 31:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod31;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod31;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss31;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss31;
                break;
            case 41:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod41;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod41;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss41;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss41;
                break;
            case 51:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod51;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod51;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss51;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss51;
                break;
            case 61:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod61;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod61;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss61;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss61;
                break;
            case 71:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod71;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod71;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss71;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss71;
                break;
            case 81:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod81;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod81;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss81;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss81;
                break;
            case 91:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod91;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod91;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss91;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss91;
                break;
            case 101:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod101;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod101;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss101;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss101;
                break;
            case 121:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod121;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod121;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss121;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss121;
                break;
            case 201:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod201;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod201;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss201;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss201;
                break;
            default:
                gaussKronrodAbscissae = QuadratureKronrodValuesType::abscissaeGaussKronrod15;
                kronrodWeights = QuadratureKronrodValuesType::weightsGaussKronrod15;
                gaussAbscissae = QuadratureKronrodValuesType::abscissaeGauss15;
                gaussWeights = QuadratureKronrodValuesType::weightsGauss15;
                break;
        }
        size_t kronrodSize = gaussKronrodAbscissae.rows();
        size_t gaussSize = gaussAbscissae.rows();

        fout << "\n//Nodes and Weights - Rule " << kronrodRule << "\n";

        //Abscissae Gauss Kronrod
        fout << "template <typename Scalar>\n"
             << "Array<Scalar, " << kronrodSize << ", 1> QuadratureKronrod<Scalar>::" 
                << gaussKronrodAbscissaeNames[i]<< " =\n"
             << "  (Array<Scalar, " << kronrodSize << ", 1>() <<\n";

        for(size_t j = 0; j < kronrodSize ; j++)
        {
            fout << "    " << gaussKronrodAbscissae(j);
            if(j !=kronrodSize - 1)
                fout << ",";
            fout << "\n";
        }

        fout << "    ).finished();\n\n";

        // Weights Gauss Kronrod
        fout << "template <typename Scalar>\n"
             << "Array<Scalar, " << kronrodSize << ", 1> QuadratureKronrod<Scalar>::" 
                << gaussKronrodWeightsNames[i] << " =\n"
             << "  (Array<Scalar, " << kronrodSize << ", 1>() <<\n";

        for(size_t j = 0; j < kronrodSize ; j++)
        {
            fout << "    " << kronrodWeights(j);
            if(j != kronrodSize - 1)
                fout << ",";
            fout << "\n";
        }

        fout << "    ).finished();\n\n";

        // Abscissae Gauss
        fout << "template <typename Scalar>\n"
             << "Array<Scalar, " << gaussSize << ", 1> QuadratureKronrod<Scalar>::" 
                << gaussAbscissaeNames[i] << " =\n"
             << "  (Array<Scalar, " << gaussSize << ", 1>() <<\n";

        for(size_t j = 0; j < gaussSize ; j++)
        {
            fout << "    " << gaussAbscissae(j);
            if(j != gaussSize - 1)
                fout << ",";
            fout << "\n";
        }

        fout << "    ).finished();\n\n";

        // Abscissae Gauss
        fout << "template <typename Scalar>\n"
             << "Array<Scalar, " << gaussSize << ", 1> QuadratureKronrod<Scalar>::" 
                << gaussWeightsNames[i] << " =\n"
             << "  (Array<Scalar, " << gaussSize << ", 1>() <<\n";

        for(size_t j = 0; j < gaussSize ; j++)
        {
            fout << "    " << gaussWeights(j);
            if(j != gaussSize - 1)
                fout << ",";
            fout << "\n";
        }

        fout << "    ).finished();\n\n";

    }

    fout << "}\n";
    fout << "#endif // EIGEN_QUADRATURE_KRONROD_H\n";

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
