#ifndef __DOP853_STEPS__
#define __DOP853_STEPS__



#include "dop853_consts.h"
#include <cmath>
#include <algorithm>
#include <utility>

namespace dop853 {

inline void step1(const size_t N,double* Y1,const double H,const double* Y,const double* K1){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*a21*K1[I];  
	}
}

inline void step2(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K2){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a31*K1[I]+a32*K2[I]);
	}
}

inline void step3(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K3){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a41*K1[I]+a43*K3[I]);
	}
}

inline void step4(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K3,const double* K4){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a51*K1[I]+a53*K3[I]+a54*K4[I]);
	}
}

inline void step5(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a61*K1[I]+a64*K4[I]+a65*K5[I]);
	}
}

inline void step6(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a71*K1[I]+a74*K4[I]+a75*K5[I]+a76*K6[I]);
	}
}



inline void step7(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6,const double* K7){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a81*K1[I]+a84*K4[I]+a85*K5[I]+a86*K6[I]+a87*K7[I]);
	}
}		



inline void step8(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6,const double* K7,const double* K8){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a91*K1[I]+a94*K4[I]+a95*K5[I]+a96*K6[I]+a97*K7[I]+a98*K8[I]);
	}
}	


inline void step9(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6,const double* K7,const double* K8,
	const double* K9){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a101*K1[I]+a104*K4[I]+a105*K5[I]+a106*K6[I]+a107*K7[I]+a108*K8[I]+a109*K9[I]);
	}
}


inline void step10(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6,const double* K7,const double* K8,
	const double* K9,const double* K10){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a111*K1[I]+a114*K4[I]+a115*K5[I]+a116*K6[I]+a117*K7[I]+a118*K8[I]+a119*K9[I]+a1110*K10[I]);
	}
}


inline void step11(const size_t N,double* Y1,const double H,const double* Y,const double* K1,
	const double* K4,const double* K5,const double* K6,const double* K7,const double* K8,
	const double* K9,const double* K10,const double* K2){
	for(size_t I=0;I<N;++I){
		Y1[I]=Y[I]+H*(a121*K1[I]+a124*K4[I]+a125*K5[I]+a126*K6[I]+a127*K7[I]+a128*K8[I]+a129*K9[I]+a1210*K10[I]+a1211*K2[I]);
	}
}



inline void stepe(const size_t N,double* K4,double* K5,const double H,const double* Y,
	const double* K1,const double* K6,const double* K7,const double* K8,const double* K9,
	const double* K10,const double* K2,const double* K3){
	for(size_t I=0;I<N;++I){
		K4[I]=b1*K1[I]+b6*K6[I]+b7*K7[I]+b8*K8[I]+b9*K9[I]+b10*K10[I]+b11*K2[I]+b12*K3[I];
		K5[I]=Y[I]+H*K4[I];
	}
}


inline std::pair<double,double> error_est(const size_t N,const double RTOLI,const double ATOLI,const double* Y,
	const double* K1,const double* K2,const double* K3,const double* K4,const double* K5,
	const double* K6,const double* K7,const double* K8,const double* K9,const double* K10){

	double ERR=0.0;
	double ERR2=0.0;
	for(size_t I=0;I<N;++I){
        double SK=ATOLI+RTOLI*std::max(std::abs(Y[I]),std::abs(K5[I]));
        double ERRI=(K4[I]-bhh1*K1[I]-bhh2*K9[I]-bhh3*K3[I])/SK;
        ERR2+=ERRI*ERRI;
        ERRI=(er1*K1[I]+er6*K6[I]+er7*K7[I]+er8*K8[I]+er9*K9[I]+er10*K10[I]+er11*K2[I]+er12*K3[I])/SK;
		ERR+=ERRI*ERRI;
	}
	return std::make_pair(ERR,ERR2);
}




inline void copy(const size_t N,double* A,const double* B){
	for(int I=0;I<N;++I)
		A[I] = B[I];
	}

}
#endif