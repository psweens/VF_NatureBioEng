//
//  viscor.cpp
//  Vascular Flows 2.0
//
//  Computation of segment viscosity based on Pries et al. Circ Res. 75,
//  904-915, 1994. Diameter corrected for differences in MCV
//
//  Created by Paul Sweeney on 06/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <cmath>

#include "global_params.h"

double viscor(const double &d, const double &hd) {
    
    double dOff = 2.4;
    double dCrit = 10.5;
    double d50 = 100;
    double eAmp = 1.1;
    double eWidth = 0.03;
    double ePeak = 0.6;
    double eHD = 1.18;
    double wMax = 2.6;
    //vplas *= 1e-9; // Viscosity of plasma in kg/mu s
    
    double wAs = 0.;
    if (dOff < d)   {
        wAs = wMax*(d-dOff)/(d+d50-2*dOff);
    }
    
    double wPeak = 0.;
    if (d > dOff && d <= dCrit) {
        wPeak = eAmp*(d-dOff)/(dCrit-dOff);
    }
    else if (dCrit < d) {
        wPeak = eAmp*exp(-eWidth*(d-dCrit));
    }
    
    double wPH = wAs + wPeak*ePeak;
    double wEFF = wAs + wPeak*(1 + hd*eHD);
    double dPH = d - 2*wPH;
    
    float hdref = 0.45F;
    double C = (cpar(0) + exp(cpar(1)*dPH)) * (-1. + 1./(1. + powf(10.F,cpar(2)) * powf(dPH,cpar(3)))) + 1./(1. + powf(10.F,cpar(2)) * powf(dPH,cpar(3))); // Curvature of relationship between relative apparent viscosity and hematrocrit with tube diameter.
    double eta45 = viscpar(0) * exp(viscpar(1)*dPH) + viscpar(2) + viscpar(3) * exp(viscpar(4) * pow(dPH,viscpar(5))); // mu_0.45
    double hdfac = (powf(1.F - hd,C) - 1.)/(powf(1.F - hdref,C) - 1.); // Discharge hematocrit fraction in mu_vitro equation
    double etaVitro = 1. + (eta45 - 1.) * hdfac;
    
    double etaVivo = etaVitro*pow(d/(d-2*wEFF),4)*vplas;
    
    return etaVivo;
}