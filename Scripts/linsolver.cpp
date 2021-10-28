//
//  linsolver.cpp
//  Vascular Flows 2.0
//
//  Solving Ax = B
//  Notation is consistent with Fry et al. (2012)
//  Solve network with known boundary conditions ... lucky you
//
//  Created by Paul Sweeney on 08/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>
#include <iostream>
#include "global_params.h"

using namespace arma;
using namespace std;

void linsolver(const int &nseg, const int &nnod, const int &nnodbc, const uvec &ista, const uvec &iend, const uvec &nodtyp, uvec &bctyp, const uvec &bcnod, vec &cond, vec &bcprfl, vec &nodpress, vec &q)    {
    
    if (phaseseparation == 0 && continuum == 0)   {
        cout<<"\t\t\tConstructing system..."<<endl;
    }
    
    
    // Construct conductance, M
    M = zeros<sp_mat>(nseg,nnod);
    for (int iseg = 0; iseg < nseg; iseg++) {
        double t_cond = cond(iseg);
        M(iseg,ista(iseg)) = t_cond;
        M(iseg,iend(iseg)) = -t_cond;
    }
    
    
    
    // Construct K matrix
    K = L * M;
    
    
    // Define qo - pressure and flow boundary conditions
    Qo = zeros<vec>(nnod);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 0)  {
            Qo(bcnod(inodbc)) = bcprfl(inodbc)*alpha;
            K.row(bcnod(inodbc)).zeros();
            K(bcnod(inodbc),bcnod(inodbc)) = 1.;
        }
        else    {
            Qo(bcnod(inodbc)) = -bcprfl(inodbc)*Gamma;
        }
    }
    
    
    nodpress = spsolve(K,Qo)/alpha;
    q = (M * nodpress)*(alpha/Gamma);
}
