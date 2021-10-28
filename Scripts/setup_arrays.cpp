//
//  setupArrays.cpp
//  Vascular Flows 2.0
//
//  Initialisation of network arrays.
//
//  Created by Paul Sweeney on 03/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>

#include "global_variables.h"
#include "global_params.h"

void setup_arrays()  {
    
    // Vectors
    ista = zeros<uvec>(nseg);
    iend = zeros<uvec>(nseg);
    nk = zeros<ivec>(nseg);
    
    nodtyp = zeros<uvec>(nnod);
    nodout = zeros<ivec>(nnod);
    nodrank = zeros<ivec>(nnod);
    
    bcnod = zeros<uvec>(nnodbc);
    
    conduc = zeros<vec>(nseg);
    c = zeros<vec>(nseg);
    rseg = zeros<vec>(nseg);
    qold = zeros<vec>(nseg);
    hdold = zeros<vec>(nseg);
    qq = zeros<vec>(nseg);
    segpress = zeros<vec>(nseg);
    tau = zeros<vec>(nseg);
    vel = zeros<vec>(nseg);
    
    ss = zeros<vec>(3);
    
    // Matrices
    nodnod = zeros<imat>(nodsegm,nnod);
    nodseg = zeros<imat>(nodsegm,nnod);
    
    End = zeros<mat>(3,nseg);
    scos = zeros<mat>(3,nseg);
    start = zeros<mat>(3,nseg);
    
}
