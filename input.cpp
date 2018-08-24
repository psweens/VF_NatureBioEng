//
//  input.cpp
//  Vascular Flows 2.0
//
//  Script reads in network data file
//
//  Dimensions:
//  Pressures - mmHg
//  Flows - um^3/s
//  Viscosity - cP
//  Diameter - um
//
//  Created by Paul Sweeney on 02/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>
#include <iostream>

#include "global_variables.h"
#include "global_params.h"
#include "initial_setup.h"
#include "misc_func.h"

using namespace std;
using namespace arma;

void input()    {
    
    int max=200;
    char bb[200];

    
    FILE *ifp;
    
    string networkFile = "Network.dat";
    networkFile = loadroot + networkFile;
    
    cout<<"\nImporting network data..."<<endl;
    ifp = fopen(networkFile.c_str(),"r");
    
    fgets(bb,max,ifp);
    printf("\t %s\n",bb);
    
    netName = bb;
    
    fscanf(ifp, "%f %f %f", &alx,&aly,&alz); fgets(bb,max,ifp);
    fscanf(ifp, "%i %i %i", &mxx,&myy,&mzz); fgets(bb,max,ifp);
    fscanf(ifp, "%f", &lb); fgets(bb,max,ifp);
    fscanf(ifp, "%f", &maxl); fgets(bb,max,ifp);
    fscanf(ifp,"%i", &nodsegm);
    fgets(bb,max,ifp);
    // Number of segments in vessel network
    fscanf(ifp,"%i", &nseg); fgets(bb,max,ifp);
    fgets(bb,max,ifp);
    // Segment properties: name type nodename(start), nodename(end), diameter, flow, hematocrit
    segname = zeros<uvec>(nseg);
    vesstyp = zeros<uvec>(nseg);
    segnodname = zeros<umat>(2,nseg);
    diam = zeros<vec>(nseg);
    lseg = zeros<vec>(nseg);
    q = zeros<vec>(nseg);
    hd = zeros<vec>(nseg);
    compLseg = 0;
    int num = detect_col(ifp);
    if (num == 7)   {
        for(int iseg = 0; iseg < nseg; iseg++){
            fscanf(ifp, "%lli %lli %lli %lli %lf %lf %lf\n",
                   &segname(iseg),&vesstyp(iseg),&segnodname(0,iseg),&segnodname(1,iseg),&diam(iseg),&q(iseg),&hd(iseg));
        }
        compLseg = 1;
    }
    else if (num == 8)  {
        for(int iseg = 0; iseg < nseg; iseg++){
            fscanf(ifp, "%lli %lli %lli %lli %lf %lf %lf %lf\n",
                   &segname(iseg),&vesstyp(iseg),&segnodname(0,iseg),&segnodname(1,iseg),&diam(iseg),&lseg(iseg),&q(iseg),&hd(iseg));
        }
    }
    else    {
        printf("*** Error in Network File: Invalid Segment Format ***");
    }

    
    // Number of nodes in vessel network
    fscanf(ifp,"%i", &nnod);
    fgets(bb,max,ifp);
    fgets(bb,max,ifp);
    // Coordinates of nodes
    nodname = zeros<uvec>(nnod);
    cnode = zeros<mat>(3,nnod);
    nodpress = zeros<vec>(nnod);
    num = detect_col(ifp);
    if (num == 4)   {
        for(int inod = 0; inod < nnod; inod++)  {
            fscanf(ifp, "%lli %lf %lf %lf\n", &nodname(inod),&cnode(0,inod),&cnode(1,inod),&cnode(2,inod));
        }
    }
    else if (num == 5)  {
        for(int inod = 0; inod < nnod; inod++)  {
            fscanf(ifp, "%lli %lf %lf %lf %lf\n", &nodname(inod),&cnode(0,inod),&cnode(1,inod),&cnode(2,inod),&nodpress(inod));
        }
    }
    else    {
        printf("*** Error in Network File: Invalid Node Format ***");
    }
    
    // Boundary nodes
    fscanf(ifp,"%i", &nnodbc);
    fgets(bb,max,ifp);
    fgets(bb,max,ifp);
    
    bcnodname = zeros<uvec>(nnodbc);
    bctyp = zeros<uvec>(nnodbc);
    bcprfl = zeros<vec>(nnodbc);
    bchd = zeros<vec>(nnodbc);
    PO2 = zeros<vec>(nnodbc);
    BCflow = zeros<vec>(nnodbc);
    BCpress = zeros<vec>(nnodbc);
    
    num = detect_col(ifp);
    if (num == 4)   {
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            fscanf(ifp,"%lli %lli %lf %lf\n", &bcnodname(inodbc),&bctyp(inodbc),&bcprfl(inodbc),&bchd(inodbc));
        }
    }
    else if (num == 5)  {
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            fscanf(ifp,"%lli %lli %lf %lf %lf\n", &bcnodname(inodbc),&bctyp(inodbc),&bcprfl(inodbc),&bchd(inodbc),&PO2(inodbc));
        }
    }
    else if (num == 6)  {
        for(int inodbc = 0; inodbc < nnodbc; inodbc++){
            fscanf(ifp,"%lli %lli %lf %lf %lf %lf\n", &bcnodname(inodbc),&bctyp(inodbc),&BCpress(inodbc),&BCflow(inodbc),&bchd(inodbc),&PO2(inodbc));
        }
        bcprfl(find(bctyp == 0)) = BCpress(find(bctyp == 0));
        bcprfl(find(bctyp == 1)) = BCflow(find(bctyp == 1));
    }
    else    {
        printf("*** Error in Network File: Invalid Boundary Node Format ***");
    }
    
    fclose(ifp);
    
    nodsegm += 1; // Armadillo indexing starts at zero
    
    
    rheolParams();

    
}
