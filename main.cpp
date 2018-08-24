//
//  main.cpp
//  Vascular Flows 2.0
//
//  Flow and Pressure Solver in Discrete Vascular Networks making use of the Armadillo
//  C++ library for stability and speed of simulations (requires BLAS, LAPACK and SuperLU).
//
//  Created by Paul Sweeney on 06/05/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <iostream>
#include <armadillo>
#include <cmath>

#include "initial_setup.h"
#include "sub_menus.h"
#include "output_data.h"
#include "discrete_flow.h"


using namespace arma;
using namespace std;

/******* GLOBAL VARIABLES *******/

// Misc
string root = "/Users/paulsweeney/Documents/Vascular-Flow-4.0/Build_Data/";    // Folder to output all created files
// Requires folders 'Discrete Flow' & 'Network Geometry' in the above root
string loadroot = "/Users/paulsweeney/Documents/Vascular-Flow-4.0/";    // Address of files to be loaded
int continuum=0,compLseg,geoAlgo=0,sec_min;
double picfac,netVol,vasDen,RBF,inflow,run_time,run_start,high,low,perc1,perc2,vol,rho;
ivec dataFitFlag;
vec vel;

// 'Network.dat' file variables
string netName;

int mxx,myy,mzz,nodsegm,nseg,nnod,nnodbc;
float alx,aly,alz,lb,maxl,totalq;

uvec segname,vesstyp,nodname,bcnodname,bctyp,BCgeo;
vec diam,q,hd,bcprfl,bchd,PO2;

umat segnodname;
mat cnode,bcp;

// 'RheolParams.dat' variables
int varyviscosity,phaseseparation;

float qtol,hdtol,constvisc,vplas,mcv,consthd;

vec bifpar,cpar,viscpar;

// Array Setup variables
float tlength;

uvec ista,iend,nodtyp,bcnod;
ivec nk, nodout;
vec conduc,c,lseg,rseg,qold,hdold,qq,segpress,tau,nodpress,ss;

imat nodnod,nodseg;
mat End,scos,start;

// Variable hematocrit variable
int nnodfl;

ivec nodrank;

// Flow estimation variables
int unknod,nIBnod;
float Gamma,alpha,beta,xi;
double targPress,targStress,ktau,oldktau,kp;
ivec unknod_vec,zeroFlowVec;
vec B,Qo,p0,tau0,storeBC,storeBChd,flowSign,oldFlowSign,oldTau,BCpress,BCflow;
uvec storeBCtyp;
sp_mat A,K,L,M,H,W,Mt;


int main(int argc, char** argv) {
    
    // Check directories
    //check_dir();
    
    // Import network file
    input();

    
    // Setup network arrays
    setup_arrays();
    
    
    // Analyse imported network
    analyse_net();

    
    // Output data for convex hull 
    hull_print("Print_Hull",0);
    hull_print("Print_Hull",1);
    
    
    // Statistical analysis of network
    outputNet("Print Network.txt", nseg, nnod, nnodbc, segname, vesstyp, segnodname, diam, q, hd, nodname, cnode, bcnodname, bctyp, bcprfl, bchd);
    output_analysis("Data_Analysis.txt");
    netStatistics("Network_Statistics.txt", "Diameter", nseg, diam, "Length", nseg, lseg, 0, "Node Type", nnod, nodtyp);
    

    
    // Discrete flow menu
    cout<<"\n\t...Discrete solver options:"<<endl;
    flow_menu();
    

    
    // Interstitial Flow
    printf("\n\t...Solving interstitial pressure\n");
    inter_menu();
    
    
    printf("\n...Fin\n");
    
    return 0;
}
