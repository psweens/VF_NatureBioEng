//
//  discrete_flow.h
//  Vascular Flows 2.0
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#ifndef Vascular_Flows_2_0_discrete_flow_h
#define Vascular_Flows_2_0_discrete_flow_h


#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

extern void assignBCs(const int &numBC);
extern double viscor(const double &d, const double &h);
extern void flow(const int &nseg, const int &nnod, const int &nnodbc, const uvec &segname, const uvec &nodname, const uvec &bcnodname, const mat &cnode, const uvec &ista, const uvec &iend, vec &diam, vec &lseg, vec &q, vec &qold, vec &hd, vec &hdold, uvec &nodtyp, uvec &bctyp, const uvec &bcnod, vec &bchd, vec &bcprfl, vec &c, vec &cond, vec &segpress, vec &nodpress, vec &qq, vec &tau, vec &BCpress, vec &BCflow);
extern void linsolver(const int &nseg, const int &nnod, const int &nnodbc, const uvec &ista, const uvec &iend, const uvec &nodtyp, uvec &bctyp, const uvec &bcnod, vec &cond, vec &bcprfl, vec &nodpress, vec &q);
extern void flowest(const int &nseg, const int &nnod, const int &nnodbc, const uvec &ista, const uvec &iend, const uvec &nodtyp, uvec &bctyp, const uvec &bcnod, const vec &lseg, vec &cond, vec &c, vec &bcprfl, vec &nodpress, vec &q);
extern void flowSetup();
extern void flowArrays();
extern void flowOutput();
extern void outputFlow(const string &filename, const uvec &vesstyp, const int &flag, const string &param1, const string &param2, const vec &press, const vec &qq, const int &logScale);
extern cube discretePerfusion(const int &disc_x, const int &disc_y, const int &disc_z);
extern void randomBC(const int &check);
extern void flow_analysis(const string &filename);
extern void flow_script(const int &check);
extern void vasc_perf(const int &tag);


#endif
