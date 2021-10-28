//
//  setupArrays2.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "inter_var.h"

void setupArrays2() {
    
    
    iK = zeros<vec>(ns);
    r0 = zeros<vec>(ns);
    spress = zeros<vec>(ns);
    s_q = zeros<vec>(ns);
    
    
    xsl0 = zeros<vec>(3);
    xsl1 = zeros<vec>(3);
    xsl2 = zeros<vec>(3);
    
    
}
