//
//  global_params.h
//  Vascular Flows
//
//  Created by Paul Sweeney on 06/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#ifndef global_params_h
#define global_params_h

#include <armadillo>

using namespace arma;


// Misc
extern int continuum,compLseg,geoAlgo,sec_min;
extern double picfac,netVol,vasDen,RBF,inflow,run_time,run_start,high,low,perc1,perc2,vol,rho;


// Network file variables
extern int mxx,myy,mzz,nodsegm;
extern float alx,aly,alz,lb,maxl,totalq;


// Rheology
extern int varyviscosity,phaseseparation;
extern float qtol,hdtol,constvisc,vplas,mcv,consthd;
extern vec bifpar,cpar,viscpar;


// Variable Hd
extern int nnodfl;


// Flow estimation
extern int unknod,nIBnod;
extern float Gamma,alpha,beta,xi;
extern double targPress,targStress,ktau,oldktau,kp;
extern ivec unknod_vec,zeroFlowVec;
extern vec B,Qo,p0,tau0,storeBC,storeBChd,oldTau,BCpress,BCflow,flowSign,oldFlowSign;
extern uvec storeBCtyp;
extern sp_mat A,K,L,M,H,W,Mt;


#endif /* global_params_h */
