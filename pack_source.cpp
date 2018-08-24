//
//  circle_packing.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 26/11/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_variables.h"
#include "global_params.h"
#include "inter_var.h"
#include "misc_func.h"

void pack_source()   {
    
    r0 *= 1e3;
    snode *= 1e3;
    
    ivec flag_circle = sphere_packing(0.95, r0, snode);

    int tag = 0;
    for (int is = 0; is < ns; is++) {
        if (flag_circle(is) == 2)    {
            snode.shed_col(is);
            source_flag.shed_row(is);
            flag_circle.shed_row(is);
            s_lseg.shed_row(is);
            s_rseg.shed_row(is);
            s_segpress.shed_row(is);
            sname.shed_row(is);
            iK.shed_row(is);
            r0.shed_row(is);
            spress.shed_row(is);
            s_q.shed_row(is);
            old_q.shed_row(is);
            ns -= 1;
            is = 0;
            tag = 1;
        }
    }
    if (tag == 1)   {
        outputf(ift,"No. of Interstitial Sources (Update)", ns,"");
    }
    
    
    outputf(ift,"Min. Source Radius", min(r0),"um");
    outputf(ift,"Max. Source Radius", max(r0),"um");
    
    
    FILE *ipf;
    
    string filename = "Interstitial/r0_xyz_radii.txt";
    filename = root + filename;
    
    ipf = fopen(filename.c_str(),"w");
    
    for (int is = 0; is < ns; is++) {
        fprintf(ipf,"%f %f %f %f\n",snode(0,is),snode(1,is),snode(2,is),r0(is));
    }
    
    fclose(ipf);
    
    r0 *= 1e-3;
    snode *= 1e-3;
    
}
