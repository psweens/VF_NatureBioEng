//
//  inter_menu.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include "inter_func.h"
#include "global_variables.h"
#include "global_params.h"
#include "misc_func.h"
#include "initial_setup.h"
#include "output_data.h"

FILE *ift;

int ns;
double kappa,iP,Lp,oncotic,PI_b,PI_v,sig;
uvec source_flag,sname,tort_net;
vec r0,s_q,s_p,spress,iK,xsl0,xsl1,xsl2,s_lseg,s_rseg,s_segpress,old_q;
mat snode,G;
fmat float_G;

void inter_menu()   {
    
    for (int iseg = 0; iseg < nseg; iseg++) {
        segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
    }
    
    ivec flag = zeros<ivec>(nnodbc);
    mat perf_xyz;
    string name = "Read_Hull";
    hull_read(name,flag,perf_xyz,1);
    
    int load_inter = 0;
    if (load_inter == 0)    {
        
        // Initialise Greens log
        init_log(ift);
        
        run_start = clock();
        
        // Locate sources
        find_sources();
        
        // Flow solver
        iflow();
        
        time_check(ift, run_start, "Interstitial flow solver run-time: ");
        
    }
        
    
    // Contour plots
    DContour("Interstitial/Pressure.ps",0);
    DContour("Interstitial/Velocity.ps",1);
        
    
    // Calculate discretised perfusion field & IFV
    finite_mesh(30,29,33);
    
    
    // Calculate tissue perfusion
    tiss_perf();
    
    
        

}
