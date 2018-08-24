//
//  inter_var.h
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#ifndef inter_var_h
#define inter_var_h

#include <armadillo>

using namespace arma;

extern FILE *ift;
extern int ns;
extern double kappa,iP,Lp,oncotic,PI_b,PI_v,sig;
extern uvec source_flag,sname,tort_net;
extern vec r0,spress,s_q,s_p,iK,xsl0,xsl1,xsl2,s_lseg,s_rseg,s_segpress,old_q;
extern mat snode,G;
extern fmat float_G;


#endif /* inter_var_h */
