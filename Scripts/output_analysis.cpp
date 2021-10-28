//
//  output_analysis.cpp
//  Vascular-Flow
//
//  General network statistical analysis
//
//  Created by Paul Sweeney on 01/07/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_variables.h"
#include "global_params.h"

void output_analysis(const string &filename)  {
    
    // Network volume
    netVol = sum(lseg % (diam % diam))*M_PI/4;
    
    // Vascular density
    vasDen = 100*netVol/(alx*aly*alz);

    // Vascular Length Density
    double lden = sum(lseg)/(alx*aly*alz)*1e6; //mm^-2
    
    // Vascular surface density
    double sden = sum(diam % lseg)*M_PI/(alx*aly*alz)*1e3;
    
    // Vascular surface area to volume ratio
    double SV = sum(diam % lseg)*M_PI/netVol*1e3;
    
    // Maximum extravascular diffusion distance
    double R = 1/sqrt(M_PI*lden*1e-6);
    
    netVol /= 1.e9; // mm^3
    
    cout<<"\t Network Volume = "<<netVol<<"mm^3"<<endl;
    cout<<"\t Vascular Density = "<<vasDen<<"%"<<endl;
    cout<<"\t Max Diameter = "<<max(diam)<<"um"<<endl;
    cout<<"\t Min Diameter = "<<min(diam)<<"um"<<endl;
    cout<<"\t Max Length = "<<max(lseg)<<"um"<<endl;
    cout<<"\t Min Length = "<<min(lseg)<<"um"<<endl;
    
    FILE *ofp;
    
    string rootname = root + filename;
    
    ofp = fopen(rootname.c_str(),"w");
    fprintf(ofp,"Network Vol. = %.3f mm^3\n",netVol);
    fprintf(ofp,"Vascular Density = %.3f %%\n",vasDen);
    fprintf(ofp,"Vascular Length Density = %.3f mm^-2\n",lden);
    fprintf(ofp,"Vascular Surface Density = %.3f mm^-1\n",sden);
    fprintf(ofp,"Vascular Surface Area/Vascular Volume = %.3f mm^-1\n",SV);
    fprintf(ofp,"Max. Extravascular Diffusion Distance = %.3f um^-1\n",R);
    fprintf(ofp,"Mean Diam. = %.3f ± %.3f um\n",mean(diam),stddev(diam));
    fprintf(ofp,"Min. Diam. = %.3f um\n",diam.min());
    fprintf(ofp,"Max. Diam. = %.3f um\n",diam.max());
    fprintf(ofp,"Mean Length = %.3f ± %.3f um\n",mean(lseg),stddev(lseg));
    fprintf(ofp,"Min. Length = %.3f um\n",lseg.min());
    fprintf(ofp,"Max. Length = %.3f um\n",lseg.max());
    
    fclose(ofp);
    
    
}
