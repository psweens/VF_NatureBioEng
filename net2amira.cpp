//
//  net2amira.cpp
//  Vascular Flows 2.0
//
//  Script to output inputted network to an Amira file format
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <iostream>

#include "output_data.h"

using namespace std;

extern string root;

void net2amira(const string &filename, const string &param_name, const int &nnod, const int &nseg, const mat &cnode, const uvec &ista, const uvec &iend, const vec &rseg, const vec &param)   {
    
    FILE *ofp1;
    
    string rootname = root + filename;
    ofp1 = fopen(rootname.c_str(),"w");
    fprintf(ofp1,"# AmiraMesh 3D ASCII 2.0 \n \n ");
    fprintf(ofp1,"\n define VERTEX %i",nnod);
    fprintf(ofp1,"\n define EDGE %i",nseg);
    fprintf(ofp1,"\n define POINT %i",(2*nseg));
    fprintf(ofp1,"\n \n Parameters { \n \t ContentType \"HxSpatialGraph\" \n } \n");
    
    fprintf(ofp1,"\n VERTEX { float[3] VertexCoordinates } @1");
    fprintf(ofp1,"\n EDGE { int[2] EdgeConnectivity } @2");
    fprintf(ofp1,"\n EDGE { int NumEdgePoints } @3");
    fprintf(ofp1,"\n POINT { float [3] EdgePointCoordinates } @4");
    fprintf(ofp1,"\n POINT { float Radii } @5");
    fprintf(ofp1,"\n POINT { float %s } @6 \n \n",param_name.c_str());
    
    // Node coordinates
    fprintf(ofp1,"@1\n");
    for(int inod = 0; inod < nnod; inod++)  {
        fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,inod),cnode(1,inod),cnode(2,inod));
    }
    
    // Connecting nodes
    fprintf(ofp1,"\n@2\n");
    
    for(int iseg = 0; iseg < nseg; iseg++)  {
        fprintf(ofp1,"%i %i\n",(int) ista(iseg),(int) iend(iseg));
    }
    
    // Number of points per edge
    fprintf(ofp1,"\n@3\n");
    for(int iseg = 0; iseg < nseg; iseg++)  {
        fprintf(ofp1,"%i \n",2);
    }
    
    // Coordinates of points - in order of start/end nodes for a segment (as read in network data file)
    fprintf(ofp1,"\n@4\n");
    for(int iseg = 0; iseg < nseg; iseg++)  {
        fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,ista(iseg)),cnode(1,ista(iseg)),cnode(2,ista(iseg)));
        fprintf(ofp1,"%.15e %.15e %.15e\n",cnode(0,iend(iseg)),cnode(1,iend(iseg)),cnode(2,iend(iseg)));
    }
    
    // Segment radii
    fprintf(ofp1,"\n@5\n");
    for(int iseg = 0; iseg < nseg; iseg++)    {
        fprintf(ofp1,"%.15e\n",rseg(iseg));
        fprintf(ofp1,"%.15e\n",rseg(iseg));
    }
    
    // Segment parameter
    fprintf(ofp1,"\n@6\n");
    for(int iseg = 0; iseg < nseg; iseg++)    {
        fprintf(ofp1,"%.15e\n",param(iseg));    // Need to have as many values as points (hence two fprintf)
        fprintf(ofp1,"%.15e\n",param(iseg));
    }
    
    fclose(ofp1);
}