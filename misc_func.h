//
//  misc_func.h
//  Vascular-Flow
//
//  Created by Paul Sweeney on 04/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#ifndef misc_func_h
#define misc_func_h

#include <armadillo>

using namespace arma;

extern int detect_col(FILE *ifp);
extern void outputf(FILE *ofp, const string &text, const double &num, const string &unit);
extern void init_log(FILE *ofp);
extern int max_double(const double &a, const double &b);
extern ivec sphere_packing(const double grow_decay, vec &r0, const mat &snode);
extern void time_check(FILE *ift, double &run_start, const string &text);


#endif /* misc_func_h */
