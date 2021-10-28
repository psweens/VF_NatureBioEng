//
//  rheolParams.cpp
//  Vascular Flows
//
//  Blood rheology parameters
//
//  Created by Paul on 31/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>

#include "global_params.h"

void rheolParams()  {
    
    
    cout<<"\tAssigning rheology parameters..."<<endl;
    
    bifpar = zeros<vec>(3);
    cpar = zeros<vec>(4);
    viscpar = zeros<vec>(6);

    
    // Bifurcation law
    bifpar(0) = 0.964;
    bifpar(1) = 6.98;
    bifpar(2) = -13.29;
    
    // C - viscosity variable
    cpar(0) = 0.80;
    cpar(1) = -0.075;
    cpar(2) = -11.0;
    cpar(3) = 12.0;
    
    
    // In vitro visocity params
    viscpar(0) = 220;
    viscpar(1) = -1.3;
    viscpar(2) = 3.2;
    viscpar(3) = -2.44;
    viscpar(4) = -0.06;
    viscpar(5) = 0.645;
    
    // Constant visc
    constvisc = 3.;
    
    // Plasma viscosity (cP)
    vplas = 1.0466;
    
    // Constant hematocrit
    consthd = 0.45;
    
    // Mean cell vol.
    mcv = 55.;
    
    // Variable viscosity & phase separation - (1) yes, (0) no
    varyviscosity = 1;
    phaseseparation = 0;
    
    
    // Hd and q tolerances for variable Hd
    hdtol = 1.e-3;
    qtol = 1.e-3;
    

    
}
