//
//  OutputData.cpp
//  Vascular Flows 2.0
//
//  Outputs raw data
//
//  Created by Paul Sweeney on 17/02/2016.
//  Copyright (c) 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

extern string root;


void outputData(const string &filename, const vec &vec1, const vec &vec2)    {
    
    FILE *ofp;
    
    string rootname = root + filename;
    
    ofp = fopen(rootname.c_str(),"w");
    for (int i = 0; i < vec1.n_elem; i++)   {
        fprintf(ofp,"%f %f\n",vec1(i),vec2(i));
    }
    
    fclose(ofp);
    
}