//
//  global_variables.h
//  Vascular Flows 2.0
//
//  Created by Paul Sweeney on 02/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#ifndef __Vascular_Flows_2_0__global_variables__
#define __Vascular_Flows_2_0__global_variables__

#include <stdio.h>
#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;


// Misc
extern string root,loadroot;
extern ivec dataFitFlag;
extern vec vel;


// 'Network.dat' variables
extern string netName;

extern int nseg,nnod,nnodbc;

extern uvec segname,vesstyp,nodname,bcnodname,bctyp,BCgeo;
extern vec diam,q,hd,bcprfl,bchd,PO2;

extern umat segnodname;
extern mat cnode,bcp;


// Array Setup variables
extern float tlength;

extern uvec ista,iend,nodtyp,bcnod;
extern ivec nk, nodout;
extern vec conduc,c,lseg,rseg,qold,hdold,qq,segpress,tau,nodpress,ss;

extern imat nodnod,nodseg;
extern mat End,scos,start;


// Variable hematocrit variable
extern ivec nodrank;


#endif /* defined(__Vascular_Flows_2_0__global_variables__) */
