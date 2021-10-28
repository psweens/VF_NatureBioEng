//
//  flowest.cpp
//  Vascular Flows 2.0
//
//  Constructing flow estimation matrices in the case of unknown boundary conditions
//  Solving the system Ax = B
//  Notation is consistent with Fry et al. (2012)
//  Note, designed to reduce run-time and memory for large networks and so split into
//  sections for phase separation or constant network haematocrit
//
//  Created by Paul Sweeney on 08/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "discrete_flow.h"
#include "global_params.h"


void flowest(const int &nseg, const int &nnod, const int &nnodbc, const uvec &ista, const uvec &iend, const uvec &nodtyp, uvec &bctyp, const uvec &bcnod, const vec &lseg, vec &cond, vec &c, vec &bcprfl, vec &nodpress, vec &q)   {
    
    if (phaseseparation == 1)   {
        
        // Construct conductance, M
        M = zeros<sp_mat>(nseg,nnod);
        for (int iseg = 0; iseg < nseg; iseg++) {
            double t_cond = cond(iseg);
            M(iseg,ista(iseg)) = t_cond;
            M(iseg,iend(iseg)) = -t_cond;
        }
        
        
        
        // Construct K matrix but inputting into matrix A to reduce run-time
        int size = nnod + nIBnod;
        A = zeros<sp_mat>(size,size);
        A(span(0,nIBnod-1),span(0,nnod-1)) = L * M;
        int cntr = 0;
        for (int inod = 0; inod < nnod; inod++) {
            if (unknod_vec(inod) == 0)  {
                for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                    if (bcnod(inodbc) == inod && bctyp(inodbc) == 0)    {
                        A.row(inod - cntr).zeros();
                        A(inod - cntr,inod) = 1.;
                    }
                }
            }
            else {
                cntr += 1;
            }
        }
        
        
        
        // Constructing H matrix - replacing with matrix A for reduced run-time
        if (nseg < 1e3) {
            A(span(nIBnod,size-1),span(0,nnod-1)) = ((repmat((square(c).t()),nnod,1) % (M.t()))*(repmat(lseg,1,nnod) % M))*ktau + W;
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
            
            A(span(nIBnod,size-1),span(0,nnod-1)) = ktau*(H1 * H2) + W;
        }
        
        
        
        // Final A matrix entry - transpose of K
        A(span(nIBnod,size-1),span(nnod,size-1)) = A(span(0,nIBnod-1),span(0,nnod-1)).t();
        
        // Construct B vector
        B(span(nIBnod,size-1)) = (W*p0)+(ktau*(M.t()*(tau0 % (c % lseg))));
        
        
    }
    else    {
        
        int size = nnod + nIBnod;
        A(span(nIBnod,size-1),span(0,nnod-1)) = ktau*H + W;
        
        // Construct B vector
        B(span(nIBnod,size-1)) = (W*p0)+(ktau*(Mt*(tau0 % (c % lseg))));
        
        
    }
    
    
    
    if (phaseseparation == 0 && continuum == 0)   {
        cout<<"\t\t\tSolving system"<<endl;
    }
    
    vec x;
    if (phaseseparation == 1)   {
        
        // Equilibrates system - scales rows and columns to unit norm in order to prevent singular system
        superlu_opts settings;
        settings.equilibrate = true;
        
        x = spsolve(A,B,"superlu",settings);
        
    }
    else {
        x = spsolve(A,B);
    }
    
    
    nodpress = x(span(0,nnod-1))/alpha;
    q = (M * nodpress)*alpha/Gamma;
    
}
