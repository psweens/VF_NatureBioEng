//
//  flowArrays.cpp
//  Vascular Flows
//
//  Created by Paul Sweeney on 08/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include "global_variables.h"
#include "global_params.h"

void flowArrays()   {
    
    flowSign = zeros<vec>(nseg);
    oldFlowSign = zeros<vec>(nseg);
    tau0 = zeros<vec>(nseg);
    oldTau = zeros<vec>(nseg);
    
    BCpress = zeros<vec>(nnodbc);
    BCflow = zeros<vec>(nnodbc);
    
}
