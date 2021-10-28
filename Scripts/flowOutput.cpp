//
//  flowOutput.cpp
//  Vascular Flows
//
//  Outputs flow data files - stats, Amira etc
//
//  Created by Paul Sweeney on 23/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>

#include "global_variables.h"
#include "global_params.h"
#include "discrete_flow.h"
#include "output_data.h"


void flowOutput()   {
    
    FILE *ofp;
    
    string rootname = "Conditions.txt";
    rootname = loadroot + rootname;
    ofp = fopen(rootname.c_str(),"w");
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        fprintf(ofp,"%lf %lf %lf %i %lf\n",cnode(0,bcnod(inodbc)),cnode(1,bcnod(inodbc)),cnode(2,bcnod(inodbc)),0,BCflow(inodbc));
    }
    
    fclose(ofp);
    
    
    outputNet("Discrete_Flow/BCpress_Net_Output.txt", nseg, nnod, nnodbc, segname, vesstyp, segnodname, diam, q, hd, nodname, cnode, bcnodname, zeros<uvec>(nnodbc), BCpress, bchd);
    
    flow_analysis("Discrete_Flow/Flow_Analysis.txt");

    
    net2amira("Discrete_Flow/Press2Amira.txt","Pressure", nnod, nseg, cnode, ista, iend, rseg, segpress);
    net2amira("Discrete_Flow/Flow2Amira.txt","Flow", nnod, nseg, cnode, ista, iend, rseg, qq);
    net2amira("Discrete_Flow/DirectionalFlow2Amira.txt","Flow", nnod, nseg, cnode, ista, iend, rseg, q);
    net2amira("Discrete_Flow/Velocity2Amira.txt","Flow", nnod, nseg, cnode, ista, iend, rseg, vel);
    net2amira("Discrete_Flow/Stress2Amira.txt","Stress", nnod, nseg, cnode, ista, iend, rseg, tau);
    net2amira("Discrete_Flow/Hematocrit2Amira.txt","Hd", nnod, nseg, cnode, ista, iend, rseg, hd);
    

    
}
