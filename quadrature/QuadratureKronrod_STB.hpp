/**
 * \file QuadratureKronrod.h
 * \sa R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner, QUADPACK, A Subroutine Package
 *     for Automatic Integration, Springer Verlag, 1983.
 */

#ifndef EIGEN_QUADRATURE_KRONROD_STB_HPP
#define EIGEN_QUADRATURE_KRONROD_STB_HPP

namespace Eigen
{

/**
 * \brief The abscissae and weights are given for the interval (-1,1).
 *        Because of symmetry, only the positive abscissae and their
 *        corresponding weights are given.
 *
 * \var abscissaeGaussKronrod15  The abscissae of the  15 point kronrod rule.
 * \var abscissaeGaussKronrod21  The abscissae of the  21 point kronrod rule.
 * \var abscissaeGaussKronrod31  The abscissae of the  31 point kronrod rule.
 * \var abscissaeGaussKronrod41  The abscissae of the  41 point kronrod rule.
 * \var abscissaeGaussKronrod51  The abscissae of the  51 point kronrod rule.
 * \var abscissaeGaussKronrod61  The abscissae of the  61 point kronrod rule.
 * \var abscissaeGaussKronrod71  The abscissae of the  71 point kronrod rule.
 * \var abscissaeGaussKronrod81  The abscissae of the  81 point kronrod rule.
 * \var abscissaeGaussKronrod91  The abscissae of the  91 point kronrod rule.
 * \var abscissaeGaussKronrod101 The abscissae of the 101 point kronrod rule.
 * \var abscissaeGaussKronrod201 The abscissae of the 201 point kronrod rule.
 *
 * \var weightsGaussKronrod15  The weights of the  15 point kronrod rule.
 * \var weightsGaussKronrod21  The weights of the  21 point kronrod rule.
 * \var weightsGaussKronrod31  The weights of the  31 point kronrod rule.
 * \var weightsGaussKronrod41  The weights of the  41 point kronrod rule.
 * \var weightsGaussKronrod51  The weights of the  51 point kronrod rule.
 * \var weightsGaussKronrod61  The weights of the  61 point kronrod rule.
 * \var weightsGaussKronrod71  The weights of the  71 point kronrod rule.
 * \var weightsGaussKronrod81  The weights of the  81 point kronrod rule.
 * \var weightsGaussKronrod91  The weights of the  91 point kronrod rule.
 * \var weightsGaussKronrod101 The weights of the 101 point kronrod rule.
 * \var weightsGaussKronrod201 The weights of the 201 point kronrod rule.
 *
 * \var weightsGauss15  The weights of the   7 point gauss rule.
 * \var weightsGauss21  The weights of the  10 point gauss rule.
 * \var weightsGauss31  The weights of the  15 point gauss rule.
 * \var weightsGauss41  The weights of the  20 point gauss rule.
 * \var weightsGauss51  The weights of the  25 point gauss rule.
 * \var weightsGauss61  The weights of the  30 point gauss rule.
 * \var weightsGauss71  The weights of the  35 point gauss rule.
 * \var weightsGauss81  The weights of the  40 point gauss rule.
 * \var weightsGauss91  The weights of the  45 point gauss rule.
 * \var weightsGauss101 The weights of the  50 point gauss rule.
 * \var weightsGauss201 The weights of the 100 point gauss rule.
 */
template <typename Scalar>
class QuadratureKronrod
{
public:

    typedef Kronrod::LaurieGautschi<Scalar> LaurieGautschiPolicy;
    typedef typename LaurieGautschiPolicy::VectorType VectorType;

    /*
    QuadratureKronrod()
    {
        ComputeNodesAndWeights();
    }
    */

    static void ComputeNodesAndWeights()
    {
        if(compute)
        {

            Array<Scalar, Dynamic, 2> xwGK = Kronrod::multiPrecisionKronrod<Scalar>(7);
            Array<Scalar, Dynamic, 2> xwG = Kronrod::multiPrecisionGauss<Scalar>(7);
            for(size_t i=0;i<8;++i)
            {
                abscissaeGaussKronrod15(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod15(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<4;++i)
            {
                weightsGauss15(i) = fabs( xwG(i,1) );
            }


            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(10);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(10);
            for(size_t i=0;i<11;++i)
            {
                abscissaeGaussKronrod21(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod21(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<5;++i)
            {
                weightsGauss21(i) = fabs( xwG(i,1) );
            }


            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(15);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(15);
            for(size_t i=0;i<16;++i)
            {
                abscissaeGaussKronrod31(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod31(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<8;++i)
            {
                weightsGauss31(i) = fabs( xwG(i,1) );
            }


            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(20);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(20);
            for(size_t i=0;i<21;++i)
            {
                abscissaeGaussKronrod41(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod41(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<10;++i)
            {
                weightsGauss41(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(25);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(25);
            for(size_t i=0;i<26;++i)
            {
                abscissaeGaussKronrod51(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod51(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<13;++i)
            {
                weightsGauss51(i) = fabs( xwG(i,1) );
            }


            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(30);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(30);
            for(size_t i=0;i<31;++i)
            {
                abscissaeGaussKronrod61(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod61(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<15;++i)
            {
                weightsGauss61(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(35);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(35);
            for(size_t i=0;i<36;++i)
            {
                abscissaeGaussKronrod71(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod71(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<18;++i)
            {
                weightsGauss71(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(40);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(40);
            for(size_t i=0;i<41;++i)
            {
                abscissaeGaussKronrod81(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod81(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<20;++i)
            {
                weightsGauss81(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(45);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(45);
            for(size_t i=0;i<46;++i)
            {
                abscissaeGaussKronrod91(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod91(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<23;++i)
            {
                weightsGauss91(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(50);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(50);
            for(size_t i=0;i<51;++i)
            {
                abscissaeGaussKronrod101(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod101(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<25;++i)
            {
                weightsGauss101(i) = fabs( xwG(i,1) );
            }

            xwGK = Kronrod::multiPrecisionKronrod<Scalar>(100);
            xwG = Kronrod::multiPrecisionGauss<Scalar>(100);
            for(size_t i=0;i<101;++i)
            {
                abscissaeGaussKronrod201(i) = fabs( xwGK(i,0) );
                weightsGaussKronrod201(i) =  fabs( xwGK(i,1) );
            }
            for(size_t i=0;i<50;++i)
            {
                weightsGauss201(i) = fabs( xwG(i,1) );
            }


            //std::cout<<"computing for precision ="<<Scalar::get_default_prec()<<std::endl;
/*
            VectorType x=VectorType::Zero(15);
            VectorType w=VectorType::Zero(15);
            LaurieGautschiPolicy::mpkonrad(7,x,w);
            for(size_t i=0;i<8;++i)
            {
                abscissaeGaussKronrod15(i) = fabs( x(i) );
                weightsGaussKronrod15(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(7);
            w=VectorType::Zero(7);
            LaurieGautschiPolicy::mpgauss(7,x,w);
            for(size_t i=0;i<4;++i)
            {
                weightsGauss15(i) = fabs( w(i) );
            }



            x=VectorType::Zero(21);
            w=VectorType::Zero(21);
            LaurieGautschiPolicy::mpkonrad(10,x,w);
            for(size_t i=0;i<11;++i)
            {
                abscissaeGaussKronrod21(i) = fabs( x(i) );
                weightsGaussKronrod21(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(10);
            w=VectorType::Zero(10);
            LaurieGautschiPolicy::mpgauss(10,x,w);
            for(size_t i=0;i<5;++i)
            {
                weightsGauss21(i) = fabs( w(i) );
            }



            x=VectorType::Zero(31);
            w=VectorType::Zero(31);
            LaurieGautschiPolicy::mpkonrad(15,x,w);
            for(size_t i=0;i<16;++i)
            {
                abscissaeGaussKronrod31(i) = fabs( x(i) );
                weightsGaussKronrod31(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(15);
            w=VectorType::Zero(15);
            LaurieGautschiPolicy::mpgauss(15,x,w);
            for(size_t i=0;i<8;++i)
            {
                weightsGauss31(i) = fabs( w(i) );
            }



            x=VectorType::Zero(41);
            w=VectorType::Zero(41);
            LaurieGautschiPolicy::mpkonrad(20,x,w);
            for(size_t i=0;i<21;++i)
            {
                abscissaeGaussKronrod41(i) = fabs( x(i) );
                weightsGaussKronrod41(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(20);
            w=VectorType::Zero(20);
            LaurieGautschiPolicy::mpgauss(20,x,w);
            for(size_t i=0;i<10;++i)
            {
                weightsGauss41(i) = fabs( w(i) );
            }



            x=VectorType::Zero(51);
            w=VectorType::Zero(51);
            LaurieGautschiPolicy::mpkonrad(25,x,w);
            for(size_t i=0;i<26;++i)
            {
                abscissaeGaussKronrod51(i) = fabs( x(i) );
                weightsGaussKronrod51(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(25);
            w=VectorType::Zero(25);
            LaurieGautschiPolicy::mpgauss(25,x,w);
            for(size_t i=0;i<13;++i)
            {
                weightsGauss51(i) = fabs( w(i) );
            }



            x=VectorType::Zero(61);
            w=VectorType::Zero(61);
            LaurieGautschiPolicy::mpkonrad(30,x,w);
            for(size_t i=0;i<31;++i)
            {
                abscissaeGaussKronrod61(i) = fabs( x(i) );
                weightsGaussKronrod61(i) =  fabs( w(i) );
            }
            x=VectorType::Zero(30);
            w=VectorType::Zero(30);
            LaurieGautschiPolicy::mpgauss(30,x,w);
            for(size_t i=0;i<15;++i)
            {
                weightsGauss61(i) = fabs( w(i) );
            }
*/
            //////////////////// PRINT VALUES ////////////////////////////////
            /*
            std::cout<<"\nGaussKronrod15 \n"<<std::endl;
            std::cout<<std::fixed;
            for(size_t i=0;i<8;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod15(i)
                    <<"\t"<<weightsGaussKronrod15(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<4;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss15(i)<<std::endl;
            }

            std::cout<<"\nGaussKronrod21 \n"<<std::endl;
            for(size_t i=0;i<11;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod21(i)
                    <<"\t"<<weightsGaussKronrod21(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<5;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss21(i)<<std::endl;
            }


            std::cout<<"\nGaussKronrod31 \n"<<std::endl;
            for(size_t i=0;i<16;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod31(i)
                    <<"\t"<<weightsGaussKronrod31(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<8;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss31(i)<<std::endl;
            }


            std::cout<<"\nGaussKronrod41 \n"<<std::endl;
            for(size_t i=0;i<21;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod41(i)
                    <<"\t"<<weightsGaussKronrod41(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<10;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss41(i)<<std::endl;
            }

            std::cout<<"\nGaussKronrod51 \n"<<std::endl;
            for(size_t i=0;i<26;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod51(i)
                    <<"\t"<<weightsGaussKronrod51(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<13;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss51(i)<<std::endl;
            }


            std::cout<<"\nGaussKronrod61 \n"<<std::endl;
            for(size_t i=0;i<31;++i)
            {
                std::cout<<std::setprecision(33)<<abscissaeGaussKronrod61(i)
                    <<"\t"<<weightsGaussKronrod61(i)<<std::endl;
            }
            std::cout<<std::endl;
            for(size_t i=0;i<15;++i)
            {
                std::cout<<std::setprecision(33)<<weightsGauss61(i)<<std::endl;
            }
            */

            compute = false;
        }
    }

    static Array<Scalar, 8, 1> abscissaeGaussKronrod15;
    static Array<Scalar, 8, 1> weightsGaussKronrod15;
    static Array<Scalar, 4, 1> weightsGauss15;

    static Array<Scalar, 11, 1> abscissaeGaussKronrod21;
    static Array<Scalar, 11, 1> weightsGaussKronrod21;
    static Array<Scalar, 5, 1> weightsGauss21;

    static Array<Scalar, 16, 1> abscissaeGaussKronrod31;
    static Array<Scalar, 16, 1> weightsGaussKronrod31;
    static Array<Scalar, 8, 1> weightsGauss31;

    static Array<Scalar, 21, 1> abscissaeGaussKronrod41;
    static Array<Scalar, 21, 1> weightsGaussKronrod41;
    static Array<Scalar, 10, 1> weightsGauss41;

    static Array<Scalar, 26, 1> abscissaeGaussKronrod51;
    static Array<Scalar, 26, 1> weightsGaussKronrod51;
    static Array<Scalar, 13, 1> weightsGauss51;

    static Array<Scalar, 31, 1> abscissaeGaussKronrod61;
    static Array<Scalar, 31, 1> weightsGaussKronrod61;
    static Array<Scalar, 15, 1> weightsGauss61;

    static Array<Scalar, 36, 1> abscissaeGaussKronrod71;
    static Array<Scalar, 36, 1> weightsGaussKronrod71;
    static Array<Scalar, 18, 1> weightsGauss71;

    static Array<Scalar, 41, 1> abscissaeGaussKronrod81;
    static Array<Scalar, 41, 1> weightsGaussKronrod81;
    static Array<Scalar, 20, 1> weightsGauss81;

    static Array<Scalar, 46, 1> abscissaeGaussKronrod91;
    static Array<Scalar, 46, 1> weightsGaussKronrod91;
    static Array<Scalar, 23, 1> weightsGauss91;

    static Array<Scalar, 51, 1> abscissaeGaussKronrod101;
    static Array<Scalar, 51, 1> weightsGaussKronrod101;
    static Array<Scalar, 25, 1> weightsGauss101;

    static Array<Scalar, 101, 1> abscissaeGaussKronrod201;
    static Array<Scalar, 101, 1> weightsGaussKronrod201;
    static Array<Scalar, 50, 1> weightsGauss201;

    static bool compute;
};

template <typename Scalar>
bool QuadratureKronrod<Scalar>::compute = true;

template <typename Scalar>
Array<Scalar, 8, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod15 =
  (Array<Scalar, 8, 1>() <<
   0.991455371120812639206854697526329,
   0.949107912342758524526189684047851,
   0.864864423359769072789712788640926,
   0.741531185599394439863864773280788,
   0.586087235467691130294144838258730,
   0.405845151377397166906606412076961,
   0.207784955007898467600689403773245,
   0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 8, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod15 =
  (Array<Scalar, 8, 1>() <<
     0.022935322010529224963732008058970,
     0.063092092629978553290700663189204,
     0.104790010322250183839876322541518,
     0.140653259715525918745189590510238,
     0.169004726639267902826583426598550,
     0.190350578064785409913256402421014,
     0.204432940075298892414161999234649,
     0.209482141084727828012999174891714
  ).finished();

template <typename Scalar>
Array<Scalar, 4, 1> QuadratureKronrod<Scalar>::weightsGauss15 =
  (Array<Scalar, 4, 1>() <<
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
  ).finished();

template <typename Scalar>
Array<Scalar, 11, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod21 =
  (Array<Scalar, 11, 1>() <<
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 11, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod21 =
  (Array<Scalar, 11, 1>() <<
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821
  ).finished();

template <typename Scalar>
Array<Scalar, 5, 1> QuadratureKronrod<Scalar>::weightsGauss21 =
  (Array<Scalar, 5, 1>() <<
    0.066671344308688137593568809893332,
    0.149451349150580593145776339657697,
    0.219086362515982043995534934228163,
    0.269266719309996355091226921569469,
    0.295524224714752870173892994651338
  ).finished();

template <typename Scalar>
Array<Scalar, 16, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod31 =
  (Array<Scalar, 16, 1>() <<
    0.998002298693397060285172840152271,
    0.987992518020485428489565718586613,
    0.967739075679139134257347978784337,
    0.937273392400705904307758947710209,
    0.897264532344081900882509656454496,
    0.848206583410427216200648320774217,
    0.790418501442465932967649294817947,
    0.724417731360170047416186054613938,
    0.650996741297416970533735895313275,
    0.570972172608538847537226737253911,
    0.485081863640239680693655740232351,
    0.394151347077563369897207370981045,
    0.299180007153168812166780024266389,
    0.201194093997434522300628303394596,
    0.101142066918717499027074231447392,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 16, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod31 =
  (Array<Scalar, 16, 1>() <<
    0.005377479872923348987792051430128,
    0.015007947329316122538374763075807,
    0.025460847326715320186874001019653,
    0.035346360791375846222037948478360,
    0.044589751324764876608227299373280,
    0.053481524690928087265343147239430,
    0.062009567800670640285139230960803,
    0.069854121318728258709520077099147,
    0.076849680757720378894432777482659,
    0.083080502823133021038289247286104,
    0.088564443056211770647275443693774,
    0.093126598170825321225486872747346,
    0.096642726983623678505179907627589,
    0.099173598721791959332393173484603,
    0.100769845523875595044946662617570,
    0.101330007014791549017374792767493
  ).finished();

template <typename Scalar>
Array<Scalar, 8, 1> QuadratureKronrod<Scalar>::weightsGauss31 =
  (Array<Scalar, 8, 1>() <<
    0.030753241996117268354628393577204,
    0.070366047488108124709267416450667,
    0.107159220467171935011869546685869,
    0.139570677926154314447804794511028,
    0.166269205816993933553200860481209,
    0.186161000015562211026800561866423,
    0.198431485327111576456118326443839,
    0.202578241925561272880620199967519
  ).finished();

template <typename Scalar>
Array<Scalar, 21, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod41 =
  (Array<Scalar, 21, 1>() <<
    0.998859031588277663838315576545863,
    0.993128599185094924786122388471320,
    0.981507877450250259193342994720217,
    0.963971927277913791267666131197277,
    0.940822633831754753519982722212443,
    0.912234428251325905867752441203298,
    0.878276811252281976077442995113078,
    0.839116971822218823394529061701521,
    0.795041428837551198350638833272788,
    0.746331906460150792614305070355642,
    0.693237656334751384805490711845932,
    0.636053680726515025452836696226286,
    0.575140446819710315342946036586425,
    0.510867001950827098004364050955251,
    0.443593175238725103199992213492640,
    0.373706088715419560672548177024927,
    0.301627868114913004320555356858592,
    0.227785851141645078080496195368575,
    0.152605465240922675505220241022678,
    0.076526521133497333754640409398838,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 21, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod41 =
  (Array<Scalar, 21, 1>() <<
    0.003073583718520531501218293246031,
    0.008600269855642942198661787950102,
    0.014626169256971252983787960308868,
    0.020388373461266523598010231432755,
    0.025882133604951158834505067096153,
    0.031287306777032798958543119323801,
    0.036600169758200798030557240707211,
    0.041668873327973686263788305936895,
    0.046434821867497674720231880926108,
    0.050944573923728691932707670050345,
    0.055195105348285994744832372419777,
    0.059111400880639572374967220648594,
    0.062653237554781168025870122174255,
    0.065834597133618422111563556969398,
    0.068648672928521619345623411885368,
    0.071054423553444068305790361723210,
    0.073030690332786667495189417658913,
    0.074582875400499188986581418362488,
    0.075704497684556674659542775376617,
    0.076377867672080736705502835038061,
    0.076600711917999656445049901530102
  ).finished();

template <typename Scalar>
Array<Scalar, 10, 1> QuadratureKronrod<Scalar>::weightsGauss41 =
  (Array<Scalar, 10, 1>()  <<
    0.017614007139152118311861962351853,
    0.040601429800386941331039952274932,
    0.062672048334109063569506535187042,
    0.083276741576704748724758143222046,
    0.101930119817240435036750135480350,
    0.118194531961518417312377377711382,
    0.131688638449176626898494499748163,
    0.142096109318382051329298325067165,
    0.149172986472603746787828737001969,
    0.152753387130725850698084331955098
  ).finished();

template <typename Scalar>
Array<Scalar, 26, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod51 =
  (Array<Scalar, 26, 1>() <<
    0.999262104992609834193457486540341,
    0.995556969790498097908784946893902,
    0.988035794534077247637331014577406,
    0.976663921459517511498315386479594,
    0.961614986425842512418130033660167,
    0.942974571228974339414011169658471,
    0.920747115281701561746346084546331,
    0.894991997878275368851042006782805,
    0.865847065293275595448996969588340,
    0.833442628760834001421021108693570,
    0.797873797998500059410410904994307,
    0.759259263037357630577282865204361,
    0.717766406813084388186654079773298,
    0.673566368473468364485120633247622,
    0.626810099010317412788122681624518,
    0.577662930241222967723689841612654,
    0.526325284334719182599623778158010,
    0.473002731445714960522182115009192,
    0.417885382193037748851814394594572,
    0.361172305809387837735821730127641,
    0.303089538931107830167478909980339,
    0.243866883720988432045190362797452,
    0.183718939421048892015969888759528,
    0.122864692610710396387359818808037,
    0.061544483005685078886546392366797,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 26, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod51 =
  (Array<Scalar, 26, 1>() <<
    0.001987383892330315926507851882843,
    0.005561932135356713758040236901066,
    0.009473973386174151607207710523655,
    0.013236229195571674813656405846976,
    0.016847817709128298231516667536336,
    0.020435371145882835456568292235939,
    0.024009945606953216220092489164881,
    0.027475317587851737802948455517811,
    0.030792300167387488891109020215229,
    0.034002130274329337836748795229551,
    0.037116271483415543560330625367620,
    0.040083825504032382074839284467076,
    0.042872845020170049476895792439495,
    0.045502913049921788909870584752660,
    0.047982537138836713906392255756915,
    0.050277679080715671963325259433440,
    0.052362885806407475864366712137873,
    0.054251129888545490144543370459876,
    0.055950811220412317308240686382747,
    0.057437116361567832853582693939506,
    0.058689680022394207961974175856788,
    0.059720340324174059979099291932562,
    0.060539455376045862945360267517565,
    0.061128509717053048305859030416293,
    0.061471189871425316661544131965264,
    0.061580818067832935078759824240066
  ).finished();

template <typename Scalar>
Array<Scalar, 13, 1> QuadratureKronrod<Scalar>::weightsGauss51 =
  (Array<Scalar, 13, 1>() <<
    0.011393798501026287947902964113235,
    0.026354986615032137261901815295299,
    0.040939156701306312655623487711646,
    0.054904695975835191925936891540473,
    0.068038333812356917207187185656708,
    0.080140700335001018013234959669111,
    0.091028261982963649811497220702892,
    0.100535949067050644202206890392686,
    0.108519624474263653116093957050117,
    0.114858259145711648339325545869556,
    0.119455763535784772228178126512901,
    0.122242442990310041688959518945852,
    0.123176053726715451203902873079050
  ).finished();

template <typename Scalar>
Array<Scalar, 31, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod61 =
  (Array<Scalar, 31, 1>() <<
    0.999484410050490637571325895705811,
    0.996893484074649540271630050918695,
    0.991630996870404594858628366109486,
    0.983668123279747209970032581605663,
    0.973116322501126268374693868423707,
    0.960021864968307512216871025581798,
    0.944374444748559979415831324037439,
    0.926200047429274325879324277080474,
    0.905573307699907798546522558925958,
    0.882560535792052681543116462530226,
    0.857205233546061098958658510658944,
    0.829565762382768397442898119732502,
    0.799727835821839083013668942322683,
    0.767777432104826194917977340974503,
    0.733790062453226804726171131369528,
    0.697850494793315796932292388026640,
    0.660061064126626961370053668149271,
    0.620526182989242861140477556431189,
    0.579345235826361691756024932172540,
    0.536624148142019899264169793311073,
    0.492480467861778574993693061207709,
    0.447033769538089176780609900322854,
    0.400401254830394392535476211542661,
    0.352704725530878113471037207089374,
    0.304073202273625077372677107199257,
    0.254636926167889846439805129817805,
    0.204525116682309891438957671002025,
    0.153869913608583546963794672743256,
    0.102806937966737030147096751318001,
    0.051471842555317695833025213166723,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 31, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod61 =
  (Array<Scalar, 31, 1>() <<
    0.001389013698677007624551591226760,
    0.003890461127099884051267201844516,
    0.006630703915931292173319826369750,
    0.009273279659517763428441146892024,
    0.011823015253496341742232898853251,
    0.014369729507045804812451432443580,
    0.016920889189053272627572289420322,
    0.019414141193942381173408951050128,
    0.021828035821609192297167485738339,
    0.024191162078080601365686370725232,
    0.026509954882333101610601709335075,
    0.028754048765041292843978785354334,
    0.030907257562387762472884252943092,
    0.032981447057483726031814191016854,
    0.034979338028060024137499670731468,
    0.036882364651821229223911065617136,
    0.038678945624727592950348651532281,
    0.040374538951535959111995279752468,
    0.041969810215164246147147541285970,
    0.043452539701356069316831728117073,
    0.044814800133162663192355551616723,
    0.046059238271006988116271735559374,
    0.047185546569299153945261478181099,
    0.048185861757087129140779492298305,
    0.049055434555029778887528165367238,
    0.049795683427074206357811569379942,
    0.050405921402782346840893085653585,
    0.050881795898749606492297473049805,
    0.051221547849258772170656282604944,
    0.051426128537459025933862879215781,
    0.051494729429451567558340433647099
  ).finished();

template <typename Scalar>
Array<Scalar, 15, 1> QuadratureKronrod<Scalar>::weightsGauss61 =
  (Array<Scalar, 15, 1>() <<
    0.007968192496166605615465883474674,
    0.018466468311090959142302131912047,
    0.028784707883323369349719179611292,
    0.038799192569627049596801936446348,
    0.048402672830594052902938140422808,
    0.057493156217619066481721689402056,
    0.065974229882180495128128515115962,
    0.073755974737705206268243850022191,
    0.080755895229420215354694938460530,
    0.086899787201082979802387530715126,
    0.092122522237786128717632707087619,
    0.096368737174644259639468626351810,
    0.099593420586795267062780282103569,
    0.101762389748405504596428952168554,
    0.102852652893558840341285636705415
  ).finished();


template <typename Scalar>
Array<Scalar, 36, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod71 =
  (Array<Scalar, 36, 1>() <<
    0.999619298560587820419535019084310,
    0.997706569099600297260163139312095,
    0.993820293038909212258407522017840,
    0.987935764443851498035117089185486,
    0.980131657851340988019859903652610,
    0.970437616039229833215070482584770,
    0.958838696995843078921262194273922,
    0.945345148207827329538725985529975,
    0.930003753050706099225035358640641,
    0.912854261359317614464937063555764,
    0.893916305839049404824112140132683,
    0.873219125025222331523282349141385,
    0.850813544681091587042032801795713,
    0.826749899092225406834050612748558,
    0.801067213125705714703581843358378,
    0.773810252286912555267423009209887,
    0.745038975666406771644308176138581,
    0.714814501556628783264408631224447,
    0.683190418488156576777504272094656,
    0.650224364665890388675792808984559,
    0.615985710487221830539605765737640,
    0.580545344749764509934502008189690,
    0.543968351696258138206293282757722,
    0.506322773241488615024297555837327,
    0.467686183461529649022383330710758,
    0.428137541517814254187620613001480,
    0.387750696027842312609168369328125,
    0.346601554430813945876979834930238,
    0.304774001471050379619976316560777,
    0.262352941209296057970895200455581,
    0.219418258415018003189060275384260,
    0.176051061165989569974303656445060,
    0.132339270613416625611142082747192,
    0.088371343275659263600929433497549,
    0.044230407960476319024979107246482,
    0.000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 36, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod71 =
  (Array<Scalar, 36, 1>() <<
    0.0010255091107466680100695643729813,
    0.0028722600144707018817512634659694,
    0.0048980908903161470913448999535535,
    0.0068554872187842001348901767105931,
    0.0087480347678970122676293701710353,
    0.0106441267608036454852097462629852,
    0.0125521386316194285632856254932093,
    0.0144261486252936342970092860304529,
    0.0162497719998497925141143236527164,
    0.0180466511295587036618351050092727,
    0.0198246307319256843086463341682683,
    0.0215607290028207408443358229090198,
    0.0232418108954666350248316623639740,
    0.0248793898649789614231649811222091,
    0.0264787298392445215429030859931175,
    0.0280248592704803250286930840600065,
    0.0295073129404838050315984750836051,
    0.0309329856908925412458457741359957,
    0.0323057496748603256180147440938779,
    0.0336145496277949424410833145813728,
    0.0348507762898165786288363179516758,
    0.0360193210644325116028647289553999,
    0.0371234780367494765199105417986844,
    0.0381545539384517980284467051984874,
    0.0391053151646666430927322198779112,
    0.0399798348609348877927516063776062,
    0.0407813447585929379712539453324104,
    0.0415027911411049663610639783689415,
    0.0421380228974238161451798354350352,
    0.0426909348444938893679127877819585,
    0.0431650461201105957501609484475510,
    0.0435545454169731029850816654607840,
    0.0438541549245973085730802961302679,
    0.0440679834669359870768339530506186,
    0.0442000975258989694715444321923505,
    0.0442456657210562284321787960122035
  ).finished();

template <typename Scalar>
Array<Scalar, 18, 1> QuadratureKronrod<Scalar>::weightsGauss71 =
  (Array<Scalar, 18, 1>() <<
    0.0058834334204430849757538962401126,
    0.0136508283483614922664040029205164,
    0.0213229799114835808834379839662051,
    0.0288292601088942540487160397144849,
    0.0361101158634633805327169696475499,
    0.0431084223261702187823064593749082,
    0.0497693704013535298051996760849950,
    0.0560408162123701285783277471651010,
    0.0618736719660801888870141387886887,
    0.0672222852690869039643055087481486,
    0.0720447947725600646654619097852778,
    0.0763034571554420535386585378842262,
    0.0799649422423242629326620809850458,
    0.0830005937288565883799265282161770,
    0.0853866533920991252259439873911176,
    0.0871044469971835342433220316055409,
    0.0881405304302754629707388075930966,
    0.0884867949071042906382073877776157
  ).finished();

template <typename Scalar>
Array<Scalar, 41, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod81 =
  (Array<Scalar, 41, 1>() <<
    0.9997075592587000165212245421359747,
    0.9982377097105592003496227024205865,
    0.9952505734460727503656095401672364,
    0.9907262386994570064530543522213722,
    0.9847228398642500610293334148833128,
    0.9772599499837742626633702837129038,
    0.9683231268541499009037674880220805,
    0.9579168192137916558045409994527593,
    0.9460718371625000382018348083531018,
    0.9328128082786765333608521668452057,
    0.9181495430728988768290910926489990,
    0.9020988069688742967282533308684931,
    0.8846920087010897459691673156848162,
    0.8659595032122595038207818083546200,
    0.8459239855873107174207527017606315,
    0.8246122308333116631963202306660988,
    0.8020605661402521271654824805230670,
    0.7783056514265193876949715455064948,
    0.7533798034389421981719521086138979,
    0.7273182551899271032809964517549305,
    0.7001629774873299310306782124344381,
    0.6719566846141795483793545149614941,
    0.6427395243055799625372439916438246,
    0.6125538896679802379526124502306949,
    0.5814470658291300065298950120298545,
    0.5494671250951282020759313055295180,
    0.5166606073863837059773675235268673,
    0.4830758016861787129085665742448230,
    0.4487645136381637639152311335460006,
    0.4137792043716050015248797458037137,
    0.3781714354735889245684707420426704,
    0.3419940908257584730074924811791943,
    0.3053024417352467195394549972603792,
    0.2681521850072536811411843448085962,
    0.2305985218807194970043610362767275,
    0.1926975807013710997155168520651499,
    0.1545068793793944770927299658245506,
    0.1160840706752552084834512844080241,
    0.0774865883312828416911548661261720,
    0.0387724175060508219331934440246233,
    0.0000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 41, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod81 =
  (Array<Scalar, 41, 1>() <<
    0.00078786332389437149872027155489687,
    0.00220748573572677796216880922531121,
    0.00376522867934192207419437277688852,
    0.00527194271488547391100911400749610,
    0.00673181348520739996342079325197288,
    0.00819757638675148244956105325366150,
    0.00967540148401718791503549171298575,
    0.01113132166402750374938616623733449,
    0.01255438476851726603177494941601932,
    0.01396255986698061404257329276730789,
    0.01536132635910245297306719365346221,
    0.01673453247500258319616665389692845,
    0.01807386840881819058019116416727674,
    0.01938764589431774100483071215654748,
    0.02067904327352817531538698511576520,
    0.02193818733583309346140081267327422,
    0.02315893101337702414441542631657291,
    0.02434569018227335927008045325065085,
    0.02550021760313012760411536558132073,
    0.02661573749902468675858406158693142,
    0.02768762611106109151543416219022983,
    0.02871838684109212328774430242495038,
    0.02970892727777659464157767791486383,
    0.03065436089141152537823603492997426,
    0.03155122361911536248171493595674375,
    0.03240098250760594428516927456513374,
    0.03320404434125756040053598377665629,
    0.03395686283420980625135217300229064,
    0.03465693584349753394613505048683965,
    0.03530514470862184103888891975079795,
    0.03590160278362810442749426164234802,
    0.03644382653034092475806451684526801,
    0.03693016953404855460457701491128602,
    0.03736118002546921808817250104671656,
    0.03773680126309354415257237969171129,
    0.03805546377885242099072980354540710,
    0.03831632400517485967847680297981700,
    0.03851974174995072693620904485019618,
    0.03866555543914104039741924425967758,
    0.03875302937875238614021147436365714,
    0.03878210476428280538640259652566010
  ).finished();

template <typename Scalar>
Array<Scalar, 20, 1> QuadratureKronrod<Scalar>::weightsGauss81 =
  (Array<Scalar, 20, 1>() <<
    0.00452127709853319125847173287818533,
    0.01049828453115281361474217106727965,
    0.01642105838190788871286348488236393,
    0.02224584919416695726150432418420857,
    0.02793700698002340109848915750772108,
    0.03346019528254784739267818308641085,
    0.03878216797447201763997203129044616,
    0.04387090818567327199167468604171550,
    0.04869580763507223206143416044814639,
    0.05322784698393682435499647977226050,
    0.05743976909939155136661773091042599,
    0.06130624249292893916653799640839860,
    0.06480401345660103807455452956675273,
    0.06791204581523390382569010823192399,
    0.07061164739128677969548363085528683,
    0.07288658239580405906051068344251784,
    0.07472316905796826420018933626132467,
    0.07611036190062624237155807592249482,
    0.07703981816424796558830753428381025,
    0.07750594797842481126372396295832633
  ).finished();

template <typename Scalar>
Array<Scalar, 46, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod91 =
  (Array<Scalar, 46, 1>() <<
    0.9997682581981284418405050180852158,
    0.9986036451819366381565476769008205,
    0.9962364840026204293451187322006907,
    0.9926499984472037417486171205977353,
    0.9878892975583419803584513571282679,
    0.9819687150345405682393184736343415,
    0.9748747285998462628301124723762666,
    0.9666083103968946047364251608924781,
    0.9571916240132814263270907118608598,
    0.9466416909956290617847205969538371,
    0.9349627017082307284972098552538910,
    0.9221639367190003880974673609605273,
    0.9082668076834639697360506838224181,
    0.8932916717532417384646490514930573,
    0.8772515871930220121297284114898936,
    0.8601624759606642253390788705671344,
    0.8420485769104069582297548405764734,
    0.8229342205020863370357752600265020,
    0.8028389511668367746487910340869061,
    0.7817843125939062913123631880986028,
    0.7597981645204108933191171745559524,
    0.7369088489454903526237388485948920,
    0.7131412242754917418607367593196434,
    0.6885216807712005252320198258804382,
    0.6630817086022337980835232637632603,
    0.6368533944532233592712238459033832,
    0.6098660544938957089438198680823058,
    0.5821502125693531866809673344441771,
    0.5537406721593487633231289014900788,
    0.5246728204629160670911341004601624,
    0.4949796574981018370299319204102766,
    0.4646951239196350985796015023097485,
    0.4338568437417829443216644742818129,
    0.4025029438585419140779745085483465,
    0.3706693423597306189532296408558525,
    0.3383926542506021616434041000318730,
    0.3057127218662330432585853442781413,
    0.2726697697523775606087653916156452,
    0.2393018532047122417906413040485990,
    0.2056474897832637457197872254715439,
    0.1717480659497809096531343654184982,
    0.1376452059832530287565900414230656,
    0.1033783028321454046727667361505742,
    0.0689869801631441724904146141038117,
    0.0345134487517776694949153427795997,
    0.0000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 46, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod91 =
  (Array<Scalar, 46, 1>() <<
    0.00062429185718321579857663066059663,
    0.00174909651748535281285516187173086,
    0.00298417191546646641998503281623962,
    0.00417984253266509402229000496101711,
    0.00533932741808091123159062611497603,
    0.00650505471381696245904611633672299,
    0.00768283318432149209385511436140106,
    0.00884569677730051624106815234208503,
    0.00998469672095943272307954145417904,
    0.01111490008447190508121416984333782,
    0.01224167381140228811054606558412293,
    0.01335197006885853998622059924442502,
    0.01443876216592827376541528869012420,
    0.01550957989105782085098413853805036,
    0.01656805497063899757282455227144833,
    0.01760592823349179892606175575954300,
    0.01861770602562484479718322571024115,
    0.01960791451368152118747566973260205,
    0.02057913131743379223158230090375461,
    0.02152536104976848200630137341864843,
    0.02244212099311634372269863847029008,
    0.02333248915075203001289733158912301,
    0.02419844417232129113332510712740897,
    0.02503528547727422745626205933879306,
    0.02583923328997513115912223308380120,
    0.02661260570454035168063547182183052,
    0.02735705869931106756998238581749345,
    0.02806873537773934426971004212203570,
    0.02874437835191346477732940766106574,
    0.02938590583816742107250843684941804,
    0.02999482799045285292199349638404121,
    0.03056789199447173972249252053845443,
    0.03110225031086367084980881406978431,
    0.03159963354172971455344564253125186,
    0.03206152703486851892940119666539451,
    0.03248514058861336432416886222383085,
    0.03286795999547159469562119073429123,
    0.03321166534480586792126042755694344,
    0.03351780561085007806314526562438370,
    0.03378396062701886695872739048248957,
    0.03400789127506731448126888746030258,
    0.03419132258692411201484025624573950,
    0.03433593627332217900296412964336558,
    0.03443961425268577159747120385314938,
    0.03450034175251302542891776315122965,
    0.03451995999118589472369263667768431
  ).finished();

template <typename Scalar>
Array<Scalar, 23, 1> QuadratureKronrod<Scalar>::weightsGauss91 =
  (Array<Scalar, 23, 1>() <<
    0.00358266315528355893114302865935139,
    0.00832318929621824164573585312223385,
    0.01303110499158278432063108246968693,
    0.01767753525793759061709254666957709,
    0.02223984755057873239395075855216899,
    0.02669621396757766480567477879310753,
    0.03102537493451546716250793889376806,
    0.03520669220160901624769979826157509,
    0.03922023672930244756418718534392934,
    0.04304688070916497115169111308116694,
    0.04666838771837336526776847574165410,
    0.05006749923795202979913210247487432,
    0.05322801673126895194590404401931040,
    0.05613487875978647664392394037486976,
    0.05877423271884173857436151763183142,
    0.06113350083106652250188637053632557,
    0.06320144007381993774996373029066688,
    0.06496819575072343085382657035907570,
    0.06642534844984252808291471563910374,
    0.06756595416360753627091022387364863,
    0.06838457737866967453169209933431610,
    0.06887731697766132288200284829805579,
    0.06904182482923202011079855515940474
  ).finished();

template <typename Scalar>
Array<Scalar, 51, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod101 =
  (Array<Scalar, 51, 1>() <<
    0.9998119013643647189875276419197978,
    0.9988664044200710501854594449742185,
    0.9969443870188761783055252128455157,
    0.9940319694320907125851082004206947,
    0.9901650106696800094239266365815958,
    0.9853540840480058823090096256324894,
    0.9795872967607694293899342854774690,
    0.9728643851066920737133441046062521,
    0.9652016510661452004920914171912959,
    0.9566109552428079429977456441566221,
    0.9470940342449394437633588783353344,
    0.9366566189448779337808749472724966,
    0.9253135464748018862172950121859790,
    0.9130785566557918930897356427716571,
    0.8999598875642946178933512159425376,
    0.8859679795236130486375409824667536,
    0.8711191982041003104111663127679473,
    0.8554297694299460846113626439347575,
    0.8389125869672243958266945311598217,
    0.8215820708593359483562541108739395,
    0.8034568680504594110458313234433987,
    0.7845558329003992639053051963409912,
    0.7648956793723515225521339100940533,
    0.7444943022260685382605362526821942,
    0.7233727660225925157922194115651731,
    0.7015524687068222510895462578836557,
    0.6790533894067472714775057462352895,
    0.6558964656854393607816248640036798,
    0.6321050684450648764522644416531749,
    0.6077029271849502391803817963918329,
    0.5827128178178365323005934192860487,
    0.5571583045146500543155229096258016,
    0.5310648245227081395687144436535273,
    0.5044581449074642016514591318491412,
    0.4773633922330430362023884111334879,
    0.4498063349740387891471314677783758,
    0.4218141573106130970286963976729900,
    0.3934143118975651273942292538238173,
    0.3646338288616142772148509187551734,
    0.3355002454194373568369882572910717,
    0.3060421192291846090835232245823827,
    0.2762881937795319903276452785211302,
    0.2462669473981440259784123306068191,
    0.2160072368760417568472845326171013,
    0.1855385817227727727002858801589330,
    0.1548905899981459020716286209411095,
    0.1240927243591603717830011356721601,
    0.0931747015600861408544503776396004,
    0.0621665648194161690801623690624937,
    0.0310983383271888761123289896659492,
    0.0000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 51, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod101 =
  (Array<Scalar, 51, 1>() <<
    0.00050676166803489136800323005679419,
    0.00142011023816635705848924160910482,
    0.00242310374582057327186801995970143,
    0.00339459089289723743085662727982836,
    0.00433770360526373578354342360625292,
    0.00528690708431070081496631090522652,
    0.00624676165289346880527846119294363,
    0.00719585591375833177104770901157397,
    0.00812753366255938737868037294245957,
    0.00905390685313200073482984638318671,
    0.00997908051277197271520006896555491,
    0.01089303114914542139185373632048379,
    0.01179073210054285467626745948807511,
    0.01267806183701224440357173722638837,
    0.01355761601155914093045815780333978,
    0.01442332487228954690873605873699028,
    0.01527146381217306258703163983580950,
    0.01610536347829330632879586448333390,
    0.01692665787573854108211355546392361,
    0.01773116846670068736945143328713021,
    0.01851604246177747374224081714934718,
    0.01928332352180814825301889801162132,
    0.02003404583071943486093518979257579,
    0.02076512847171916224208023972166701,
    0.02147433820180993345171927343976448,
    0.02216298492660156610684349832667521,
    0.02283171756887400699163339995992034,
    0.02347818407891019300458609050334622,
    0.02410062587721830597666576014469811,
    0.02469990353732573949476020971483749,
    0.02527641257483270914647897860737446,
    0.02582833680912179520884486377406142,
    0.02635430681869013845472576260957721,
    0.02685489636290659736893714711520551,
    0.02733033189010964383305404350544173,
    0.02777922389243983336140685575895975,
    0.02820054057947867914629860187746193,
    0.02859466838340605773904731791010499,
    0.02896171979050103356757141731956686,
    0.02930066703511978233684077865154510,
    0.02961078494045995503025213539965659,
    0.02989233581703979245880464830873981,
    0.03014535377457083519407513020913017,
    0.03036913320016487318137849566585392,
    0.03056323770631003720546053583408409,
    0.03072784566941669020593437016581195,
    0.03086293412170997376026432332960616,
    0.03096809562366843795854210956298413,
    0.03104317399694570106876064569864784,
    0.03108828778240521910421818490634307,
    0.03110336664174957546715464493457709
  ).finished();

template <typename Scalar>
Array<Scalar, 25, 1> QuadratureKronrod<Scalar>::weightsGauss101 =
  (Array<Scalar, 25, 1>() <<
    0.00290862255315514095840072434285548,
    0.00675979919574540150277887817798503,
    0.01059054838365096926356968149924102,
    0.01438082276148557441937890892732435,
    0.01811556071348939035125994342235462,
    0.02178024317012479298159206906269034,
    0.02536067357001239044019487838544272,
    0.02884299358053519802990637311323243,
    0.03221372822357801664816582732300395,
    0.03545983561514615416073461100097580,
    0.03856875661258767524477015023638593,
    0.04152846309014769742241197896406702,
    0.04432750433880327549202228683039420,
    0.04695505130394843296563301363498768,
    0.04940093844946631492124358075143273,
    0.05165570306958113848990529584009528,
    0.05371062188899624652345879725566455,
    0.05555774480621251762356742561226950,
    0.05718992564772838372302931506599316,
    0.05860084981322244583512243663084847,
    0.05978505870426545750957640531258523,
    0.06073797084177021603175001538481100,
    0.06145589959031666375640678608391538,
    0.06193606742068324338408750978083069,
    0.06217661665534726232103310736061343
  ).finished();

template <typename Scalar>
Array<Scalar, 101, 1> QuadratureKronrod<Scalar>::abscissaeGaussKronrod201 =
  (Array<Scalar, 101, 1>() <<
    0.9999525032523487419455875958687268,
    0.9997137267734412336782284693423007,
    0.9992281658838012560346868946290257,
    0.9984919506395958184001633591863492,
    0.9975136131227297392523159490732886,
    0.9962951347331251491861317322411310,
    0.9948326219369267821422998436539298,
    0.9931249370374434596520098928487835,
    0.9911749876510258446843411143088782,
    0.9889843952429917480044187458077366,
    0.9865520156031485854321432946337390,
    0.9838775407060570154961001555110082,
    0.9809628425806126209479438330385846,
    0.9778093584869182885537810884292019,
    0.9744169257913283019115695267641970,
    0.9707857757637063319308978578975054,
    0.9669175365712510964078221684669607,
    0.9628136542558155272936593260301664,
    0.9584745254644430098923420791196814,
    0.9539007829254917428493369308943576,
    0.9490940500776430880700266106323840,
    0.9440558701362559779627747064152187,
    0.9387870435308086115958141514345454,
    0.9332885350430795459243336681308625,
    0.9275620589048801751310904131546791,
    0.9216092981453339526669513284819875,
    0.9154313832787401756921788026450870,
    0.9090295709825296904671263377891461,
    0.9024057102083787390268111337693575,
    0.8955616449707269866985210224302278,
    0.8884987940099135709863303098787153,
    0.8812186793850184155733168254278056,
    0.8737233059188845423290172135986970,
    0.8660146884971646234107399696762430,
    0.8580945065652815825676935927999503,
    0.8499645278795912842933625914201047,
    0.8416269226315942280342849082868531,
    0.8330838798884008235429158338447557,
    0.8243373193478551155538206639457735,
    0.8153892383391762543939887586492580,
    0.8062419752299814989885811994777252,
    0.7968978923903144763895728821832460,
    0.7873591329821532712191509680798932,
    0.7776279096494954756275513868344901,
    0.7677067275239078952210683529091099,
    0.7575981185197071760356679644384008,
    0.7473044344001021672254328341835497,
    0.7368280898020207055124277148201010,
    0.7261717523301654045878554370215473,
    0.7153381175730564464599671227043660,
    0.7043297320110017113015827506599276,
    0.6931491993558019659486479416754373,
    0.6817993432991379075019986645713245,
    0.6702830156031410158025870143232266,
    0.6586029440214330855482255981923165,
    0.6467619085141292798326303044586304,
    0.6347628808669128796170068668230022,
    0.6226088602037077716041908451723122,
    0.6103027422339497305230401468066605,
    0.5978474702471787212648065451493406,
    0.5852461548472222375133040734824663,
    0.5725019326213811913168704435257254,
    0.5596178538598728700448388327457235,
    0.5465970120650941674679942571817499,
    0.5334426463281092439192593426094116,
    0.5201580198817630566468157494552085,
    0.5067463240599214601391067274797347,
    0.4932107892081909335693087934493340,
    0.4795547716885330549048043305588901,
    0.4657816497733580422492166233957546,
    0.4518947420866361433048328724615129,
    0.4378974021720315131089780436221960,
    0.4237930916938821018150669271527927,
    0.4095852916783015425288684000571577,
    0.3952774340481667161036767355786164,
    0.3808729816246299567633625488695874,
    0.3663754887636280023186351955882702,
    0.3517885263724217209723438295489706,
    0.3371156254451936501038051146407429,
    0.3223603439005291517224765823983254,
    0.3075263156275365712075278830398810,
    0.2926171880384719647375558882354944,
    0.2776365767518349561634681429088114,
    0.2625881203715034791689293362549821,
    0.2474755186979889209661488249174269,
    0.2323024818449739696495099632079641,
    0.2170726954172689985470672220029671,
    0.2017898640957359972360488595303965,
    0.1864577395843781906516839351866999,
    0.1710800805386032748875323747070898,
    0.1556606277572599523002823353224754,
    0.1402031372361139732075146046824055,
    0.1247113982607281204528715673545208,
    0.1091892035800611150034260065793849,
    0.0936403342835439809898845668294524,
    0.0780685828134366366948173712015526,
    0.0624777614692301001035754804658728,
    0.0468716824215916316149239129338483,
    0.0312541520838667808185743504525497,
    0.0156289844215430828722166999974293,
    0.0000000000000000000000000000000000
  ).finished();

template <typename Scalar>
Array<Scalar, 101, 1> QuadratureKronrod<Scalar>::weightsGaussKronrod201 =
  (Array<Scalar, 101, 1>() <<
    0.00012796430957024721771296604719778,
    0.00035867672428027546451819698699003,
    0.00061229953852751686968468056548071,
    0.00085841483623935311664167374912298,
    0.00109796416132059890255397465155434,
    0.00133983780660382341754863019304499,
    0.00158539114055099209687612715239440,
    0.00182937001345660807993702363601705,
    0.00207023170435975873786445531287183,
    0.00231123008487149478622574422551013,
    0.00255360720035108161107860930247813,
    0.00279496754149865554132924610296587,
    0.00303416468796449379917327454133990,
    0.00327287902049999951685865323103972,
    0.00351196946843139471347132584863839,
    0.00375002082356796252973315848571211,
    0.00398619744931235184550297813723408,
    0.00422153240339578242265345682700445,
    0.00445663655824494510224150713533979,
    0.00469055184548748439481506164708745,
    0.00492264171642177778021256522028597,
    0.00515360688570317554082069990644640,
    0.00538389880957824664615043517317290,
    0.00561281343110871597561315793834040,
    0.00583984643867452187979812998241984,
    0.00606550320973303863921033440183559,
    0.00629012776973767504951202123451844,
    0.00651317411315535135623979587712922,
    0.00673423042383435906546110501890771,
    0.00695367647707533778618333444634080,
    0.00717178053315253341524750684836420,
    0.00738810306985931173687401288774170,
    0.00760230005117430829817610049263740,
    0.00781466520662034505085711719016354,
    0.00802541132258648471308962366547617,
    0.00823417501731623690079536275494453,
    0.00844066393687280440093304148919421,
    0.00864511004867877912656139959329108,
    0.00884768427998831020007201725385440,
    0.00904808029993671228835071165243691,
    0.00924604652391111724687570119177701,
    0.00944176898875185552314903500699259,
    0.00963538626179462391106485411435149,
    0.00982663642373088065571326334283804,
    0.01001530101161305053268537450183723,
    0.01020153095088495986015919422291445,
    0.01038543929706851271926198711586115,
    0.01056679983247565098373215943168121,
    0.01074542170778693080646692567379385,
    0.01092142841500597232724852424011181,
    0.01109491257094606883971232155924911,
    0.01126567746332956237923986454364446,
    0.01143355579967040516909625306002457,
    0.01159864925643444521982733110763019,
    0.01176103385541580107532667409174482,
    0.01192053786750440393445381605750602,
    0.01207701451438569994891221026442586,
    0.01223054787113266210109559371404226,
    0.01238120033058132478094426900828574,
    0.01252882177740552594442502987847792,
    0.01267328363219971069040254068582595,
    0.01281465559402401820734387481377863,
    0.01295298874891947508181321347643330,
    0.01308815204304083947131249628348671,
    0.01322003331459425763799456688719462,
    0.01334869039477997219532795119897057,
    0.01347416490042760653129384849738332,
    0.01359634288062617797837336787649308,
    0.01371512721021576932690771893822620,
    0.01383056582631592119521912896609004,
    0.01394269234050930765052410813271065,
    0.01405140838748433908894290183169631,
    0.01415663080386553141797847993528262,
    0.01425839919822706151742984801432832,
    0.01435674034589556885653893978408636,
    0.01445157028450311174994356610608770,
    0.01454281897629966555537908242478331,
    0.01463051894733967823881185958389905,
    0.01471469106719224107616070689745859,
    0.01479526485272091523273713296861773,
    0.01487218274551550851114896175149253,
    0.01494547117901571779314143269144248,
    0.01501514584992601775929409587145444,
    0.01508114903502136704084960300677971,
    0.01514343516330660012211995270464270,
    0.01520202535514477439125881646888569,
    0.01525693069928626681152511683204471,
    0.01530810567790395821753937559915117,
    0.01535551634411844871530936744265760,
    0.01539917910754267008730746329865897,
    0.01543910086714820098073247575245750,
    0.01547524789208464777691069202608381,
    0.01550759760659890369589537474066319,
    0.01553616215856523293822391298200718,
    0.01556094454341676462938645995315419,
    0.01558192251438053678351080364024222,
    0.01559908471168710670698089639871095,
    0.01561243933403481989468460959619688,
    0.01562198563724487346851127872839699,
    0.01562771265700169003321136142695118,
    0.01562962018460484993210603370275515
  ).finished();

template <typename Scalar>
Array<Scalar, 50, 1> QuadratureKronrod<Scalar>::weightsGauss201 =
  (Array<Scalar, 50, 1>() <<
    0.00073463449050567173040632065833034,
    0.00170939265351810523952935837149120,
    0.00268392537155348241943959042900112,
    0.00365596120132637518234245872752520,
    0.00462445006342211935109578908297848,
    0.00558842800386551515721194634843921,
    0.00654694845084532276415210333149526,
    0.00749907325546471157882874401639778,
    0.00844387146966897140262083490230100,
    0.00938041965369445795141823766081212,
    0.01030780257486896958578210172783538,
    0.01122511402318597711722157336633358,
    0.01213145766297949740774479244874817,
    0.01302594789297154228555858375890179,
    0.01390771070371877268795414910800464,
    0.01477588452744130176887998752035426,
    0.01562962107754600272393686595379193,
    0.01646808617614521264310498008821078,
    0.01729046056832358243934419836674167,
    0.01809594072212811666439075142049303,
    0.01888373961337490455294116588154323,
    0.01965308749443530586538147024544407,
    0.02040323264620943276683885165758377,
    0.02113344211252764154267230044096968,
    0.02184300241624738631395374130439802,
    0.02253122025633627270179697093167396,
    0.02319742318525412162248885418272729,
    0.02384096026596820596256041190228343,
    0.02446120270795705271997502334977289,
    0.02505754448157958970376422562092326,
    0.02562940291020811607564200986215087,
    0.02617621923954567634230874175730189,
    0.02669745918357096266038466418633635,
    0.02719261344657688013649156780217069,
    0.02766119822079238829420415587042646,
    0.02810275565910117331764833018699455,
    0.02851685432239509799093676286445787,
    0.02890308960112520313487622813451527,
    0.02926108411063827662011902349564095,
    0.02959048805991264251175451067883659,
    0.02989097959333283091683680666859583,
    0.03016226510516914491906868161047923,
    0.03040407952645482001650785981882518,
    0.03061618658398044849645944326205319,
    0.03079837903115259042771390303055976,
    0.03095047885049098823406346347074793,
    0.03107233742756651658781017024291803,
    0.03116383569620990678381832121718665,
    0.03122488425484935773237649864809813,
    0.03125542345386335694764247438619803
  ).finished();
}

#endif // EIGEN_QUADRATURE_KRONROD_STB_HPP
