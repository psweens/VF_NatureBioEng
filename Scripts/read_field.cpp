//
//  read_field.cpp
//  Vascular-Flow
//
//  Code to read in interstitial field data
//
//  Created by Paul Sweeney on 16/06/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "inter_var.h"
#include "global_variables.h"
#include "global_params.h"
#include "inter_func.h"

void read_field(int flag, int ziter, int nsl1, int nsl2, int NL)   {
    
    double xmin = min(cnode.row(0));
    double xmax = max(cnode.row(0));
    double ymin = min(cnode.row(1));
    double ymax = max(cnode.row(1));
    
    mat zv = zeros<mat>(nsl1+1,nsl2+1);
    vec x = zeros<vec>(3);
    
    double var_min = 1e16;
    double var_max = 1e-16;
    if (flag == 1)  {
        var_max = log(var_max);
    }
    for (int i = 0; i < ziter; i++) {
        
        string slice_name;
        slice_name = root + "Interstitial/Slices";
        if (flag == 0)  {
            slice_name += "/Pressure/Data/";
        }
        else if (flag == 1) {
            slice_name += "/Velocity/Data/";
        }
        char slice[64];
        sprintf(slice,"%i.txt",i);
        slice_name += slice;
        
        FILE *ofp = fopen(slice_name.c_str(),"r");
        float num = 0.0;
        for (int x = 0; x <= nsl1; x++) {
            for (int y = 0; y <= nsl2; y++) {
                fscanf(ofp,"%e\t",&num);
                if (num < var_min)  {
                    var_min = num;
                }
                if (num > var_max)  {
                    var_max = num;
                }
            }
            fscanf(ofp,"\n");
        }
        
        fclose(ofp);
    }
    cout<<var_min<<"\t"<<var_max<<endl;
    
    for (int i = 0; i < ziter; i++) {
        
        zv.zeros();
        string slice_name;
        slice_name = root + "Interstitial/Slices";
        if (flag == 0)  {
            slice_name += "/Pressure/Data/";
        }
        else if (flag == 1) {
            slice_name += "/Velocity/Data/";
        }
        char slice[64];
        sprintf(slice,"%i.txt",i);
        slice_name += slice;
        
        FILE *ofp = fopen(slice_name.c_str(),"r");
        for (int x = 0; x <= nsl1; x++) {
            for (int y = 0; y <= nsl2; y++) {
                fscanf(ofp,"%le\t",&zv(x,y));
            }
            fscanf(ofp,"\n");
        }
        fclose(ofp);
        
        
        cout<<i<<endl;
        char sim[64];
        sprintf(sim,"%i",i);
        string name = root + "Interstitial/Slices/";
        if (flag == 0)  {
            name += "/Pressure/";
        }
        else {
            name += "/Velocity/";
        }
        name += sim;
        name += ".ps";
        
        double pint = 0.0;
        double vmin = 0.0;
        double vmax = 0.0;
        if (flag == 1)  {
            vmin = var_min;
            vmax = var_max;
        }
        else {
            if (iP < segpress.min())    {
                vmin = iP;
            }
            else {
                vmin = segpress.min();
            }
            vmax = segpress.max();
        }
        pint = (vmax - vmin) / NL;
        cout<<vmin<<"\t"<<vmax<<endl;
        
        
        vec cl = zeros<vec>(NL+1);
        FILE *ofp2 = fopen(name.c_str(), "w");
        fprintf(ofp2, "%%!PS-Adobe-2.0\n");
        fprintf(ofp2, "/Times-Roman findfont\n");
        fprintf(ofp2, "12 scalefont\n");
        fprintf(ofp2, "setfont\n");
        for(int i = 0; i <= NL; i++) cl(i) = vmin + (i)*pint;
        
        
        DContour_shade(ofp2,nsl1+1,nsl2+1,picfac,NL,xmin,xmax,ymin,ymax,cl,zv);
        
        fclose(ofp2);
        
    }
    
    
}
