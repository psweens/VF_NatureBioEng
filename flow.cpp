//
//  flow.cpp
//  Vascular Flows 2.0
//
//  Calculates blood flow .. duh
//
//  Created by Paul Sweeney on 06/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <cmath>

#include "discrete_flow.h"
#include "output_data.h"
#include "global_params.h"

using namespace std;

extern string loadroot;

void flow(const int &nseg, const int &nnod, const int &nnodbc, const uvec &segname, const uvec &nodname, const uvec &bcnodname, const mat &cnode, const uvec &ista, const uvec &iend, vec &diam, vec &lseg, vec &q, vec &qold, vec &hd, vec &hdold, uvec &nodtyp, uvec &bctyp, const uvec &bcnod, vec &bchd, vec &bcprfl, vec &c, vec &cond, vec &segpress, vec &nodpress, vec &qq, vec &tau, vec &BCpress, vec &BCflow) {
    
    BCpress = zeros<vec>(nnodbc);
    BCflow = zeros<vec>(nnodbc);

    // Millimetre scaling
    diam *= 1e-3;
    lseg *= 1e-3;

    // Detecting unknown boundary conditions
    int fryCntr = 0;
    if (any(bctyp == 3))    {
        fryCntr = 0;
    }


    
    // Nonlinear iterations.  If convergence is slow, hematocrit is increasingly underrelaxed.
    // This eventually forces convergence even when system is unstable due to rheological
    // effects (see work of R.L. Carr et al.).
    int nitmax = 1;
    uword errsegq = 0,errseghd = 0;
    if (phaseseparation == 1)   {
        nitmax = 100;
    }
    float relax = 1., maxqerr = 0., maxhderr = 0.;
    vec qchange,hdchange;
    for (int iter = 1; iter <= nitmax; iter++)  {
        if (iter % 5 == 0)  {
            relax *= 0.8;
        }
        
        qold = q;
        hdold = hd;
        double visc = 0.;
        if (phaseseparation == 1 || accu(bctyp(find(bctyp == 3))) == 0)   {
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
                cond(iseg) = M_PI*pow(diam(iseg),4)/(128*visc*lseg(iseg));
                c(iseg) = 4*visc/(M_PI*pow((diam(iseg)*0.5),3));
            }
        }
        
        if (any(bctyp == 3))    {
            flowest(nseg, nnod, nnodbc, ista, iend, nodtyp, bctyp, bcnod, lseg, cond, c, bcprfl, nodpress, q);
        }
        else    {
            linsolver(nseg, nnod, nnodbc, ista, iend, nodtyp, bctyp, bcnod, cond, bcprfl, nodpress, q);
        }
        
        // Detecting zero flow
        int zeroFlow = 0;
        zeroFlowVec = zeros<ivec>(nseg);
        zeroFlowVec(find(q == 0.0)).ones();
        zeroFlow = (int) accu(zeroFlowVec);
        if (zeroFlow > 0 && continuum == 0)   {
            cout<<"\t\t\t*** Zero flow detected in "<<zeroFlow<<" segments ***"<<endl;
            //net2amira("Zero2Amira.txt","ZeroFlow", nnod, nseg, cnode, ista, iend, diam*1e3, conv_to<vec>::from(zeroFlowVec));
        }
        
        
        fryCntr += 1;
        
    
        
        // Compare hd and q with previous values
        maxqerr = 0.;
        maxhderr = 0.;
        errsegq = 0;
        errseghd = 0;
        qchange = q - qold;
        hdchange = hd - hdold;
        hd = hdold + relax*hdchange;
        
        maxqerr = abs(qchange).max(errsegq);
        maxhderr = abs(hdchange).max(errseghd);
        printf("\t\t\t\tIteration: %i, Flow error = %f, h'crit error = %f\n",iter,maxqerr,maxhderr);
        if(maxqerr < qtol && maxhderr < hdtol)  {
            goto converged;
        }
        
    }
    errsegq = (int) errsegq;
    errseghd = (int) errseghd;

    if (phaseseparation == 1)   {
        printf("\n\t\t\t\t*** Warning: Nonlinear iteration not converged ***\n");
        printf("\t\t\t\tFlow error = %f at segment %lli, h'crit error = %f at segment %lli\n",maxqerr,segname(errsegq),maxhderr,segname(errseghd));
    }
    converged:;
    
    // Return to micrometre scaling
    diam *= 1e3;
    lseg *= 1e3;
    
    
    // Calculate absolute flow, shear stress and segment pressure
    qq = abs(q);
    tau = (c % qq)*(Gamma/beta);
    
    
}

