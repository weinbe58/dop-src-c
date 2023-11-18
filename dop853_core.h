#ifndef __DOP853_CORE__
#define __DOP853_CORE__

#include <iostream>

#include <cmath>
#include <algorithm>
#include <utility>
#include "dop853_steps.h"
#include "dop853_consts.h"



namespace dop853 {


typedef void (*functype)(const size_t,const double,const double*,double*);

void complex_wrapper()

class dop853_core
{
    const functype FCN;
    const size_t N;
    const double ATOLI;
    const double RTOLI;
    const size_t NMAX;
    const double HMAX;
    const double SAFE;
    const double FAC1;
    const double FAC2;
    const double BETA;
    const double UROUND;

    double X;
    double H;
    double H0;
    size_t NFCN;
    bool initial_state_set;

    double* WORK;
    double* Y;
    double* Y1;
    double* K1;
    double* K2;
    double* K3;
    double* K4;
    double* K5;
    double* K6;
    double* K7;
    double* K8;
    double* K9;
    double* K10;

public:
    dop853_core(const functype,const size_t,const double,const double,const size_t,
        const double,const double,const double,const double,const double,const double);
    ~dop853_core();
    double get_X(){return X;}
    void set_initial_state(const double,const double[]);
    void get_Y(double[]);
    int integrate(const double);



};




void dop853_core::set_initial_state(const double X0,const double Y0[]){
    for(size_t I=0;I<N;++I){
        Y[I] = Y0[I];
    }
    X=X0;
    initial_state_set=true;
    NFCN=0;
}
void dop853_core::get_Y(double YOUT[]){
    for(size_t I=0;I<N;++I){
        YOUT[I] = Y[I];
    }
}


dop853_core::dop853_core(const functype _FCN,const size_t _N,
    const double _RTOLI=1.0E-6,const double _ATOLI=1.0E-12,const size_t _NMAX=100000,const double _HMAX=0.0,
    const double _H0=0.0,const double _SAFE=0.9,const double _FAC1=0.3,const double _FAC2=6.0,const double _BETA=0.0) : \
    FCN(_FCN), N(_N), RTOLI(_RTOLI), ATOLI(_ATOLI), NMAX(_NMAX), HMAX(_HMAX), H0(_H0), SAFE(_SAFE), \
    FAC1(_FAC1), FAC2(_FAC2), BETA(_BETA), UROUND(2.3E-16) {


    initial_state_set = false;
    WORK = new double[12*_N];
    Y   = WORK;
    Y1  = Y  + _N;
    K1  = Y1 + _N;
    K2  = K1 + _N;
    K3  = K2 + _N;
    K4  = K3 + _N;
    K5  = K4 + _N;
    K6  = K5 + _N;
    K7  = K6 + _N;
    K8  = K7 + _N;
    K9  = K8 + _N;
    K10 = K9 + _N;
}
dop853_core::~dop853_core(){
    delete[] WORK;
    WORK = NULL;
    Y   = NULL;
    Y1  = NULL;
    K1  = NULL;
    K2  = NULL;
    K3  = NULL;
    K4  = NULL;
    K5  = NULL;
    K6  = NULL;
    K7  = NULL;
    K8  = NULL;
    K9  = NULL;
    K10 = NULL;
}


int dop853_core::integrate(const double XEND){
// *** *** *** *** *** *** ***
//  INITIALISATIONS
// *** *** *** *** *** *** ***
    if(!initial_state_set){
        return -5;
    }
    if(X==XEND){
        return 0;
    }

    int IDID=0;
    size_t NSTEP=0;
    const double EXPO1=1.0/8.0-BETA*0.2;
    const double FACC1=1.0/FAC1;
    const double FACC2=1.0/FAC2;
    double FACOLD=1.0E-4;
    double ERR,ERR2,HNEW,XPH,FAC,FAC11,DENO;
    const double POSNEG=std::copysign(1.0,XEND-X);
//  INITIAL PREPARATIONS      
    bool LAST=false;
    bool REJECT=false;
    double HMAX_NEW;
    if(HMAX==0){
        HMAX_NEW=XEND-X;
    }
    else{
        HMAX_NEW=HMAX;
    }
    const double HMAX_ABS=std::abs(HMAX_NEW);
    if(std::abs(H0)>HMAX_ABS){
        return -6;
    }
    // step size estimate 
    double DNF=0;
    double DNY=0;
    double DER2=0;

	std::pair<double,double> ERRS;

    (*FCN)(N,X,Y,K1);
    if(H0==0){
// ----------------------------------------------------------
// ----  COMPUTATION OF AN INITIAL STEP SIZE GUESS
// ----------------------------------------------------------

// ---- COMPUTE A FIRST GUESS FOR EXPLICIT EULER AS
// ----   H = 0.01 * NORM (Y0) / NORM (F0)
// ---- THE INCREMENT FOR EXPLICIT EULER IS SMALL
// ---- COMPARED TO THE SOLUTION

        for(size_t I=0;I<N;++I){
            double SK=ATOLI+RTOLI*std::abs(Y[I]);
            DNF+=(K1[I]/SK)*(K1[I]/SK);
            DNY+=(Y[I]/SK)*(Y[I]/SK);
        }


        if(DNF<=1.0E-10 || DNY<=1.0E-10){
            H=1.0E-6;
        }
        else{
            H=std::sqrt(DNY/DNF)*0.01;
        }
        H=std::min(H,HMAX_ABS);
        H=std::copysign(H,POSNEG);

// ---- PERFORM AN EXPLICIT EULER STEP
        for(size_t I=0;I<N;++I){
            Y1[I]=Y[I]+H*K1[I];
        }

        (*FCN)(N,X+H,Y1,K2);

// ---- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION
        for(size_t I=0;I<N;++I){
            double SK=ATOLI+RTOLI*std::abs(Y[I]);
            double DIFF=(K2[I]-K1[I])/SK;
            DER2+=DIFF*DIFF;
        }


        DER2=std::sqrt(DER2)/H;
// ---- STEP SIZE IS COMPUTED SUCH THAT
// ---- H**IORD * MAX ( NORM (F0), NORM (DER2)) = 0.01
        double DER12=std::max(std::abs(DER2),std::sqrt(DNF));
        double H1;
        if(DER12<=1.0E-15){
            H1=std::max(1.0E-6,std::abs(H)*1E-3);
        }
        else{
            H1=std::pow((0.01/DER12),(1.0/8.0));
        }
        H=std::min(100.0*std::abs(H),std::min(H1,HMAX_ABS));
        H=std::copysign(H,POSNEG);
   
    }
    else{
        H=std::copysign(std::abs(H0),POSNEG);
    }

//  BASIC INTEGRATION STEP
    while(true){
        if(NSTEP>NMAX){ // IF (NSTEP.GT.NMAX) GOTO 78
            IDID = -2;
            break;
        }
        if(0.1*std::abs(H)<=std::abs(X)*UROUND){ // IF (0.1D0*ABS(H).LE.ABS(X)*UROUND)GOTO 77
            IDID = -3;
            break;
        }
        if((X+1.01*H-XEND)*POSNEG>0){
            H=XEND-X;
            LAST=true;
        }

    	NSTEP++;
// ---- THE TWELVE STAGES
        step1(N,Y1,H,Y,K1);
        (*FCN)(N,X+c2*H,Y1,K2);
        step2(N,Y1,H,Y,K1,K2);
        (*FCN)(N,X+c3*H,Y1,K3);
        step3(N,Y1,H,Y,K1,K3);
        (*FCN)(N,X+c4*H,Y1,K4);
        step4(N,Y1,H,Y,K1,K3,K4);
        (*FCN)(N,X+c5*H,Y1,K5);
        step5(N,Y1,H,Y,K1,K4,K5);
        (*FCN)(N,X+c6*H,Y1,K6);
        step6(N,Y1,H,Y,K1,K4,K5,K6);
        (*FCN)(N,X+c7*H,Y1,K7);
        step7(N,Y1,H,Y,K1,K4,K5,K6,K7);
        (*FCN)(N,X+c8*H,Y1,K8);
        step8(N,Y1,H,Y,K1,K4,K5,K6,K7,K8);
        (*FCN)(N,X+c9*H,Y1,K9);
        step9(N,Y1,H,Y,K1,K4,K5,K6,K7,K8,K9);
        (*FCN)(N,X+c10*H,Y1,K10);
        step10(N,Y1,H,Y,K1,K4,K5,K6,K7,K8,K9,K10);
        (*FCN)(N,X+c11*H,Y1,K2);
        step11(N,Y1,H,Y,K1,K4,K5,K6,K7,K8,K9,K10,K2);
        (*FCN)(N,X+H,Y1,K3);
        // NFCN+=11;
        stepe(N,K4,K5,H,Y,K1,K6,K7,K8,K9,K10,K2,K3);
// ---- ERROR ESTIMATION
        ERRS=error_est(N,RTOLI,ATOLI,Y,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10);
        ERR = ERRS.first;
        ERR2 = ERRS.second;

		DENO=ERR+0.01*ERR2;
        if(DENO<=0.0){DENO=1.0;}
        ERR=std::abs(H)*ERR*std::sqrt(1.0/(N*DENO));
// ---- COMPUTATION OF HNEW
        FAC11=std::pow(ERR,EXPO1);
// ---- LUND-STABILIZATION
        FAC=FAC11/std::pow(FACOLD,BETA);
// ---- WE REQUIRE  FAC1 <= HNEW/H <= FAC2
        FAC=std::max(FACC2,std::min(FACC1,FAC/SAFE));
        HNEW=H/FAC;

        if(ERR<=1.0){
// -------- STEP IS ACCEPTED
           	FACOLD=std::max(ERR,1.0E-4);

            (*FCN)(N,X+H,K5,K4);
// -------- STIFFNESS DETECTION (ADD LATER)
// -------- DENSE OUTPUT (NOT IMPLEMENTING)
            copy(N,K1,K4);
            copy(N,Y,K5);
// -------- NORMAL EXIT
            if(LAST){
                H=HNEW;
                break;
            }

            if(std::abs(HNEW)>HMAX_ABS){
                HNEW = POSNEG*HMAX_ABS;
            }
            if(REJECT){
                HNEW=POSNEG*std::min(std::abs(HNEW),std::abs(H));
            }
            REJECT=false;  
        	X+=H;

        }
        else
        {
// -------- STEP IS REJECTED
            HNEW=H/std::min(FACC1,FAC11/SAFE);
            LAST=false;
            REJECT=true;
        }
    	H=HNEW;
    	ERR=ERR2=0;
    }
   	H0=H; // store current step size in case user continues integration
   
    return IDID;
}



}
#endif