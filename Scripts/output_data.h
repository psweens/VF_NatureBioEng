//
//  output_data.h
//  Vascular Flows 2.0
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#ifndef Vascular_Flows_2_0_output_data_h
#define Vascular_Flows_2_0_output_data_h

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

ivec histogram(const vec &param);
extern void outputNet(const string &name, const int &nseg, const int &nnod, const int &nnodbc, const uvec &segname, const uvec &vesstyp, const umat &segnodname, const vec &diam, const vec &q, const vec &hd, const uvec &nodname, const mat &cnode, const uvec &bcnodname, const uvec &bctyp, const vec &bcprfl, const vec &bchd);
extern void netStatistics(const string &name, const string &p1name, const int &p1Total, const vec &p1, const string &p2name, const int &p2Total, const vec &p2, const int &logScale, const string &p3name, const int &p3Total, const uvec &p3);
extern void net2amira(const string &filename, const string &param_name, const int &nnod, const int &nseg, const mat &cnode, const uvec &ista, const uvec &iend, const vec &diam, const vec &param);
extern void outputData(const string &filename, const vec &vec1, const vec &vec2);
extern void output_analysis(const string &filename);
extern void hull_print(const string &filename, int typ);
extern void hull_read(const string &filename, ivec &flag, mat &coords, int typ);


// Misc
extern void resetNet(ivec &index);


#endif
