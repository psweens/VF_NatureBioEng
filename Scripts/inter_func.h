//
//  inter_menu.h
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#ifndef inter_menu_h
#define inter_menu_h

#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;

extern void setupArrays2();
extern void find_sources();
extern void iflow();
extern double eval_G(const double &r0, const double &kappa, const vec &x, const vec &y);
extern double eval_deltaG(const double &r0, double const &kappa, const vec &x, const vec &y);
extern void DContour(const string &fname, const int &flag);
extern void DContour_shade(FILE *ofp, const int &m, const int &n, double &scalefac, int &nl, const double &xmin, const double &xmax, const double &ymin, const double &ymax, const vec &cl, const mat &zv);
extern double DCeval(vec &x, const int &flag);
extern void pack_source();
extern void tiss_perf();
extern ivec tort_l();
extern void read_field(int flag, int ziter, int nsl1, int nsl2, int NL);
extern void finite_mesh(const int &disc_x, const int &disc_y, const int &disc_z);
extern double bicgstab(fmat &a, vec &b, vec &x, int &n, double &eps, int &itmax);
extern double l_bicgstab(vec &b, vec &x, int &n, double &eps, int &itmax);
extern double bicgstab_sym(vec &b, vec &x, int &n, double &eps, int &itmax);

#endif /* inter_menu_h */
