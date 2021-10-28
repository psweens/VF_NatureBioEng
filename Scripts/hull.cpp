//
//  convex_hull.cpp
//  Vascular-Flow
//
//  Functions to print and read network hull data
//
//  Created by Paul Sweeney on 01/07/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_variables.h"
#include "global_params.h"

// Prints boundary node and all node data to be read into Matlab which produces convex hull output
void hull_print(const string &filename, int typ)  {
    
    FILE *ofp;
    
    if (typ == 0)   {
        string rootname = root + filename + ".txt";
        
        ofp = fopen(rootname.c_str(),"w");
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            fprintf(ofp,"%llu %lf %lf %lf\n",bcnodname(inodbc),cnode(0,bcnod(inodbc)),cnode(1,bcnod(inodbc)),cnode(2,bcnod(inodbc)));
        }

    }
    else {
        
        string rootname = root + filename + "_Full.txt";
        
        ofp = fopen(rootname.c_str(),"w");
        for (int inod = 0; inod < nnod; inod++) {
            fprintf(ofp,"%llu %lf %lf %lf\n",nodname(inod),cnode(0,inod),cnode(1,inod),cnode(2,inod));
        }

        
    }
    
    
    fclose(ofp);
    
}

// flag - flag boundary nodes on the surface of the hull
// coords - coordinates along the hull of the vascular network
void hull_read(const string &filename, ivec &flag, mat &coords, int typ) {
    
    
    FILE *ofp;
    
    if (typ == 0)   {
        
        string rootname = loadroot + filename + ".txt";
        ofp = fopen(rootname.c_str(),"r");
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            fscanf(ofp,"%lli\n",&flag(inodbc));
        }
        
    }
    else {
        
        string rootname = loadroot + filename + "_Full.txt";
        
        ofp = fopen(rootname.c_str(),"r");
        fscanf(ofp,"%le\n",&vol);
        int num;
        fscanf(ofp,"%i\n",&num);
        int rows = num*2;
        coords = zeros<mat>(3,rows);
        for (int inod = 0; inod < rows; inod++) {
            fscanf(ofp,"%lf %lf %lf\n",&coords(0,inod),&coords(1,inod),&coords(2,inod));
        }
        
    }
    
    fclose(ofp);
    
}
