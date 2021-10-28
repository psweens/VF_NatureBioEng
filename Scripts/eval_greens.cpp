//
//  eval_greens.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>

#include "inter_var.h"

double eval_G(const double &r0, const double &kappa, const vec &x, const vec &y)   {
    
    double G = 0.0;
    
    double r = sqrt(pow(x(0) - y(0),2) + pow(x(1) - y(1),2) + pow(x(2) - y(2),2));
    
    if (r <= r0)    {
        G = (3 - pow(r/r0,2)) / (8*M_PI*kappa*r0);
    }
    else if (r > r0)    {
        G = 1 / (4*M_PI*kappa*r);
    }
    else   {
        printf("\n*** Error in G: Failure computing distance between sources ***\n");
    }
    
    return G;
}


double eval_deltaG(const double &r0, double const &kappa, const vec &x, const vec &y) {
    
    double dG = 0.0;
    
    double r = sqrt(pow(x(0) - y(0),2) + pow(x(1) - y(1),2) + pow(x(2) - y(2),2));
    
    if (r <= r0)    {
        dG = -r / (4*M_PI*kappa*pow(r0,3));
    }
    else if (r > r0)    {
        dG = -1 / (4*M_PI*kappa*pow(r,2));
    }
    else   {
        printf("\n*** Error in dG: Failure computing distance between sources ***\n");
    }
    
    return dG;
    
}
