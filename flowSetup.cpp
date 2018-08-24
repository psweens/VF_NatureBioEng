//
//  flowSetup.cpp
//  Vascular Flows
//
//  Pre-simulation set up of matrices to solve for network flow
//  The ordering of matrices is designed to reduce run-time in the case of
//  unknown boundary conditions and/or phase separation
//  Notation is consistent with Fry et al. (2012)
//
//
//  Created by Paul Sweeney on 06/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_params.h"
#include "global_variables.h"
#include "discrete_flow.h"

void flowSetup()    {
    
    cout<<"\t\t\tConstructing system..."<<endl;
    
    // Dimensional constants - for millimetre scaling
    Gamma = 1/(1e3*60);
    alpha = 0.1333;
    beta = 1e-4;
    xi = 0.001*1e-3;
    diam *= 1e-3;
    lseg *= 1e-3;
    
    kp = 0.1;
    ktau = 1e-4;
    
    
    // Target pressure
    targPress = mean(bcprfl(find(bctyp == 0))); // Set target pressure as the mean of the assign boundary pressure conditions
    cout<<"\t\t\tTarget Pressure: "<<targPress<<" mmHg"<<endl;
    p0 = zeros<vec>(nnod);
    p0.fill(alpha*targPress);
    
    
    // Target shear stress - intially set with random directions unless network flow is known
    targStress = 5.;
    int ran = 0;
    srand((double) time(0));
    for (int iseg = 0; iseg < nseg; iseg++)    {
    rewind:;
        ran = rand()%2;
        if (ran == 0)   {
            ran = -1;
        }
        tau0(iseg) = targStress*beta*ran;
    }
    if (accu(q) > 0.0)  {
        tau0 = abs(tau0);
        tau0 %= sign(q);
        cout<<"\t\t\t*** Wall shear stress initialised using network file flow direction ***"<<endl;
    }
    
    
    // Creating an index to flag boundary nodes with unknown boundary conditions
    unknod_vec = zeros<ivec>(nnod);
    unknod_vec.elem(bcnod.elem(find(bctyp == 3))).ones();
    unknod = (int) accu(unknod_vec);
    nIBnod = nnod - unknod; // No. of internal and known boundary nodes
    
    
    // Constructing vector Qo
    int cntr = 0;
    Qo = zeros<vec>(nIBnod);
    for (int inod = 0; inod < nnod; inod++) {
        if (unknod_vec(inod) == 0)  {
            for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                if (inod == bcnod(inodbc) && bctyp(inodbc) == 0) {
                    Qo(inod - cntr) = -bcprfl(inodbc)*alpha;
                }
                else if (inod == bcnod(inodbc) && (bctyp(inodbc) == 1 || bctyp(inodbc) == 2))  {
                    Qo(inod - cntr) = bcprfl(inodbc)*Gamma;
                }
            }
        }
        else {
            cntr += 1;
        }
    }
    
    
    // Construct W matrix
    W = zeros<sp_mat>(nnod,nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        long t_ista = ista(iseg);
        long t_iend = iend(iseg);
        double t_lseg = lseg(iseg);
        W(t_ista,t_ista) += (0.5*kp*t_lseg);
        W(t_iend,t_iend) += (0.5*kp*t_lseg);
    }
    
    
    // Construct L matrix
    cntr = 0;
    L = zeros<sp_mat>(nIBnod,nseg);
    for (int inod = 0; inod < nnod; inod++)    {
        if (unknod_vec(inod) == 0)  {
            int countIndx = 0;
            for (int iseg = 0; iseg < nseg; iseg++)    {
                if (inod == ista(iseg))   {
                    L(inod-cntr,iseg) = -1.;
                    countIndx += 1;
                }
                else if (inod == iend(iseg))  {
                    L(inod-cntr,iseg) = 1.;
                    countIndx += 1;
                }
                if (nodtyp(inod) == countIndx)  {
                    iseg = nseg;
                }
            }
        }
        else    {
            cntr += 1;
        }
    }
    
    
    // Partially construct B vector (for system Ax = B)
    int size = nnod + nIBnod;
    B = zeros<vec>(size);
    B(span(0,nIBnod-1)) = -Qo;
    
    
    // Subroutine for flow estimation and constant network haematocrit - run-time efficiency
    if (any(bctyp == 3) && phaseseparation == 0)    {
        
        // Construct conductance, M
        hd.fill(consthd);
        double visc = 0.;
        M = zeros<sp_mat>(nseg,nnod);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (varyviscosity == 1) {
                visc = viscor((diam(iseg)*1e3),hd(iseg))*xi;
                if (isnan(visc))    {
                    printf("\t\tViscosity error at segment %lli, h'crit = %f\n",segname(iseg),hd(iseg));
                }
            }
            else {
                visc = constvisc*xi;
            }
            c(iseg) = 4*visc/(M_PI*pow((diam(iseg)*0.5),3));
            M(iseg,ista(iseg)) = M_PI*pow(diam(iseg),4)/(128*visc*lseg(iseg));
            M(iseg,iend(iseg)) = -M_PI*pow(diam(iseg),4)/(128*visc*lseg(iseg));
        }
        Mt = M.t();
        
        
        // Construct K matrix
        K = L * M;
        cntr = 0;
        for (int inod = 0; inod < nnod; inod++) {
            if (unknod_vec(inod) == 0)  {
                for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                    if (bcnod(inodbc) == inod && bctyp(inodbc) == 0)    {
                        K.row(inod - cntr).zeros();
                        K(inod - cntr,inod) = 1.;
                    }
                }
            }
            else {
                cntr += 1;
            }
        }
        
        
        
        
        // Constructing H matrix
        if (nseg < 1e3) {
            H = ((repmat((square(c).t()),nnod,1) % (M.t()))*(repmat(lseg,1,nnod) % M));
        }
        else    {
            sp_mat H1 = zeros<sp_mat>(nnod,nseg);
            sp_mat H2 = zeros<sp_mat>(nseg,nnod);
            for (int iseg = 0; iseg < nseg; iseg++) {
                double t_lseg = lseg(iseg);
                double t_c = pow(c(iseg),2);
                long t_sta = ista(iseg);
                long t_end = iend(iseg);
                double t_Msta = M(iseg,t_sta);
                double t_Mend = M(iseg,t_end);
                H2(iseg,t_sta) = t_Msta * t_lseg;
                H1(t_sta,iseg) = t_Msta * t_c;
                H2(iseg,t_end) = t_Mend * t_lseg;
                H1(t_end,iseg) = t_Mend * t_c;
            }
            
            H = H1 * H2;
        }
        
        
        
        // Partially construct A matrix
        A = zeros<sp_mat>(size,size);
        A(span(0,nIBnod-1),span(0,nnod-1)) = K;
        A(span(nIBnod,size-1),span(nnod,size-1)) = K.t();
        
        
    }
    
    
    diam *= 1e3;
    lseg *= 1e3;
}

