//
//  outputFlow.cpp
//  Vascular Flows
//
//  Outputs raw data - handy for post-simulation processing
//
//  Created by Paul Sweeney on 23/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>
#include <iostream>

#include "output_data.h"

using namespace arma;
using namespace std;

void outputFlow(const string &filename, const uvec &vesstyp, const int &flag, const string &param1, const string &param2, const vec &press, const vec &qq, const int &logScale)   {
    
    int sum = 0;
    for (int iseg = 0; iseg < vesstyp.n_elem; iseg++) {
        if (vesstyp(iseg) == flag)  {
            sum += 1;
        }
    }

    
    string root = "Discrete_Flow/";
    
    string stats = " Stats.txt";
    string statsFile = root + filename + " " + param1 + " & " + param2 + stats;
    netStatistics(statsFile, filename + " " + param1, sum, press(find(vesstyp == flag)), filename + " " + param2, sum, qq(find(vesstyp == flag)), logScale, "Unused", 1, ones<uvec>(1));
    
    string data = " Data.txt";
    string dataFile = root + filename + " " + param1 + " & " + param2 + data;
    outputData(dataFile,press(find(vesstyp == flag)),qq(find(vesstyp == flag)));
    
    
    
}
