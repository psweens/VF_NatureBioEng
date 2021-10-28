//
//  DCeval.cpp
//  Vascular Flows 2.0
//
//  Taken from Greens O2 by Tim Secomb
//  and update for the armadillo library
//
//
//  Created by Paul Sweeney on 16/07/2015.
//

#include <stdio.h>
#include "global_variables.h"
#include "inter_var.h"

double DCeval(vec &x, const int &flag)  {
    
    double dist,pgreen = 0.0;
    
    // Initialize the pressure field to zero
    double p = 0.;
    
    // Add contributions from the arteriolar outflows
    
    x *= 1e-3;
    
    if (flag == 0)  {
        
        for (int i = 0; i < ns; i++)    {
            
            dist = sqrt(pow(x(0) - snode(0,i),2) + pow(x(1) - snode(1,i),2) + pow(x(2) - snode(2,i),2));
            
            if (dist <= r0(i))    {
                pgreen = 1*(3 - pow(dist/r0(i),2)) / (8*M_PI*kappa*r0(i));
            }
            else if (dist > r0(i))    {
                pgreen = 1 / (4*M_PI*kappa*dist);
            }
            
            p += pgreen*s_q(i); // greens function * outflows
            
        }
        
        p /= 0.1333;
        
    }
    else if (flag == 1) {
        
        for (int i = 0; i < ns; i++)    {
            
            dist = sqrt(pow(x(0) - snode(0,i),2) + pow(x(1) - snode(1,i),2) + pow(x(2) - snode(2,i),2));
            
            if (dist <= r0(i))    {
                pgreen = -dist / (4*M_PI*kappa*pow(r0(i),3));
            }
            else if (dist > r0(i))    {
                pgreen = -1 / (4*M_PI*kappa*pow(dist,2));
            }
            
            p -= kappa*pgreen*s_q(i); // greens function * outflows
            
        }
        
        // Absolute fluid velocity
        p = abs(p);  // mm/s
        
    }
    
    x *= 1e3;

    return p;
}
