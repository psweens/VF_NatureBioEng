//
//  find_sources.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 06/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>

#include "inter_func.h"
#include "inter_var.h"
#include "global_variables.h"
#include "global_params.h"
#include "misc_func.h"

void find_sources() {
    
    // Currently one source per vessel - will need to link each sources and segments if we discretise more
    ns = nseg;
    
    // Millimetre scaling
    
    // Define source locations
    snode = zeros<mat>(3,ns);
    for (int iseg = 0; iseg < nseg; iseg++) {
        snode(0,iseg) = 0.5*(cnode(0,ista(iseg)) + cnode(0,iend(iseg)));
        snode(1,iseg) = 0.5*(cnode(1,ista(iseg)) + cnode(1,iend(iseg)));
        snode(2,iseg) = 0.5*(cnode(2,ista(iseg)) + cnode(2,iend(iseg)));
    }
    
    sname = zeros<uvec>(ns);
    sname = segname;
    old_q = zeros<vec>(ns);

    
    int yup = 0;
    if (yup == 1)    {
        
        // Reduce no. of sources for tortuous networks
        s_lseg = lseg;
        s_rseg = rseg;
        s_segpress = segpress;
        source_flag = conv_to<uvec>::from(tort_l());
        tort_net = source_flag;
        uvec branch = find_unique(source_flag);
        ivec flag = zeros<ivec>(ns);
        for (int i = 0; i < branch.n_elem; i++) {
            ivec nod = zeros<ivec>(ns);
            for (int is = 0; is < ns; is++) {
                if (source_flag(is) == source_flag(branch(i)))   {
                    nod(ista(is)) += 1;
                    nod(iend(is)) += 1;
                }
            }
            
            int num = 1;
            if (num > 1)    {
                for (int is = 0; is < ns; is++) {
                    if (nod(ista(is)) == 1)   {
                        int sta = (int) iend(is);
                        flag(is) = 1;
                        int cntr = 0;
                        double cnt_lseg = 0.0;
                        for (int js = 0; js < ns; js++) {
                            if (sta == ista(js) && flag(js) == 0 && source_flag(js) == source_flag(branch(i)))    {
                                sta = (int) iend(js);
                                cntr += 1;
                                if (cntr == num)  {
                                    flag(js) = 1;
                                    cntr = 0;
                                    s_lseg(js) = cnt_lseg;
                                    cnt_lseg = 0.0;
                                }
                                else {
                                    flag(js) = 2;
                                    cnt_lseg += s_lseg(js);
                                }
                                js = 0;
                            }
                            else if (sta == iend(js) && flag(js) == 0 && source_flag(js) == source_flag(branch(i)))  {
                                sta = (int) ista(js);
                                cntr += 1;
                                if (cntr == num)  {
                                    flag(js) = 1;
                                    cntr = 0;
                                    s_lseg(js) = cnt_lseg;
                                    cnt_lseg = 0.0;
                                }
                                else {
                                    flag(js) = 2;
                                    cnt_lseg += s_lseg(js);
                                }
                                js = 0;
                            }
                        }
                    }
                    else if (nod(iend(is)) == 1)   {
                        int sta = (int) ista(is);
                        flag(is) = 1;
                        int cntr = 0;
                        double cnt_lseg = 0.0;
                        for (int js = 0; js < ns; js++) {
                            if (sta == ista(js) && flag(js) == 0 && source_flag(js) == source_flag(branch(i)))    {
                                sta = (int) iend(js);
                                cntr += 1;
                                if (cntr == num)  {
                                    flag(js) = 1;
                                    cntr = 0;
                                    s_lseg(js) = cnt_lseg;
                                    cnt_lseg = 0.0;
                                }
                                else {
                                    flag(js) = 2;
                                    cnt_lseg += s_lseg(js);
                                }
                                js = 0;
                            }
                            else if (sta == iend(js) && flag(js) == 0 && source_flag(js) == source_flag(branch(i)))  {
                                sta = (int) ista(js);
                                cntr += 1;
                                if (cntr == num)  {
                                    flag(js) = 1;
                                    cntr = 0;
                                    s_lseg(js) = cnt_lseg;
                                    cnt_lseg = 0.0;
                                }
                                else {
                                    flag(js) = 2;
                                    cnt_lseg += s_lseg(js);
                                }
                                js = 0;
                            }
                        }
                    }
                }
            }
            
        }
        
        // Remove any overlapping sources
        ivec tag = zeros<ivec>(nseg);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (tag(iseg) == 0) {
                for (int jseg = 0; jseg < nseg; jseg++)  {
                    if (iseg != jseg && tag(jseg) == 0)   {
                        if ((ista(iseg) == ista(jseg) && iend(iseg) == iend(jseg)) || (ista(iseg) == iend(jseg) && iend(iseg) == ista(jseg))  )  {
                            flag(jseg) = 2;
                            tag(jseg) = 1;
                        }
                    }
                }
            }
        }
        
        for (int is = 0; is < ns; is++) {
            if (flag(is) == 2)    {
                snode.shed_col(is);
                source_flag.shed_row(is);
                flag.shed_row(is);
                s_lseg.shed_row(is);
                s_rseg.shed_row(is);
                s_segpress.shed_row(is);
                sname.shed_row(is);
                old_q.shed_row(is);
                ns -= 1;
                is = 0;
            }
        }
        
        for (int is = 0; is < ns; is++) {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (sname(is) == segname(iseg)) {
                    old_q(is) = qq(iseg);
                }
            }
        }
        
    }
    else if (yup == 0) {
        
        // Subdivide segments if needed
        maxl = 50.;
        ivec flag = zeros<ivec>(nseg);
        s_lseg = zeros<vec>(ns);
        s_rseg = zeros<vec>(ns);
        s_segpress = zeros<vec>(ns);
        source_flag = segname;
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (lseg(iseg) > maxl)   {
                flag(iseg) = 1;
                double sub_seg = ceil(lseg(iseg) / maxl) - 1;
                source_flag.insert_rows(ns,sub_seg);
                snode.insert_cols(ns,sub_seg);
                s_lseg.insert_rows(ns,sub_seg);
                s_rseg.insert_rows(ns,sub_seg);
                s_segpress.insert_rows(ns,sub_seg);
                flag.insert_rows(ns,sub_seg);
                sname.insert_rows(ns,sub_seg);
                old_q.insert_rows(ns,sub_seg);
                double x = cnode(0,iend(iseg)) - cnode(0,ista(iseg));
                double y = cnode(1,iend(iseg)) - cnode(1,ista(iseg));
                double z = cnode(2,iend(iseg)) - cnode(2,ista(iseg));
                double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
                double psi = acos(z/r);
                double theta = atan2(y,x);
                r = 0.;
                for (int i = 0; i < sub_seg; i++)   {
                    sname(ns + i) = segname(iseg);
                    old_q(ns + i) = qq(iseg);
                    s_lseg(ns + i) = lseg(iseg) / sub_seg;
                    s_rseg(ns + i) = rseg(iseg);
                    s_segpress(ns + i) = segpress(iseg);
                    source_flag(ns + i) = source_flag(iseg);
                    r += 0.5*s_lseg(ns + i);
                    snode(0,ns + i) = r*sin(psi)*cos(theta) + cnode(0,ista(iseg));
                    snode(1,ns + i) = r*sin(psi)*sin(theta) + cnode(1,ista(iseg));
                    snode(2,ns + i) = r*cos(psi) + cnode(2,ista(iseg));
                    r += 0.5*s_lseg(ns + i);
                }
                s_lseg(iseg) = lseg(iseg) / sub_seg;
                ns += sub_seg;
            }
            else {
                s_lseg(iseg) = lseg(iseg);
            }
            
            s_rseg(iseg) = rseg(iseg);
            s_segpress(iseg) = segpress(iseg);
        }
        
        // Remove any overlapping sources
        ivec tag = zeros<ivec>(nseg);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (tag(iseg) == 0) {
                for (int jseg = 0; jseg < nseg; jseg++)  {
                    if (iseg != jseg && tag(jseg) == 0)   {
                        if ((ista(iseg) == ista(jseg) && iend(iseg) == iend(jseg)) || (ista(iseg) == iend(jseg) && iend(iseg) == ista(jseg))  )  {
                            flag(jseg) = 1;
                            tag(jseg) = 1;
                        }
                    }
                }
            }
        }
        
        int temp = ns;
        for (int is = 0; is < temp; is++) {
            if (flag(is) == 1)    {
                snode.shed_col(is);
                source_flag.shed_row(is);
                flag.shed_row(is);
                s_lseg.shed_row(is);
                s_rseg.shed_row(is);
                s_segpress.shed_row(is);
                sname.shed_row(is);
                old_q.shed_row(is);
                temp -= 1;
                ns -= 1;
                is = 0;
            }
        }
        
        
    }
    else    {
        s_lseg = lseg;
        s_rseg = rseg;
        for (int iseg = 0; iseg < nseg; iseg++) {
            segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
        }
        s_segpress = segpress;
        source_flag = conv_to<uvec>::from(tort_l());
    }
    outputf(ift,"No. of interstitial sources", (double) ns,"");
    
    
    // Scale snode to mm
    snode *= 1e-3;
    
    
    // Setup arrays
    setupArrays2();
    
    
    // Define r0 - initial r0 based on vessel radii or discretisation of vessel segment
    r0 = s_rseg*1e-3;
    // Algorithm updates r0 to prevent overlapping from neighbouring spheres
    pack_source();
    
    
    
    // Calculate intravascular resistance to fluid transport
    Lp = (2.80 / 0.1333) * 1e-6; // mm2 s2 / kg
    for (int is = 0; is < ns; is++) {
        double surf = 2*M_PI*s_rseg(is)*s_lseg(is)*1e-6;
        iK(is) = 1 / (Lp * surf);      // kg / mm2 s
    }
   
    
    outputf(ift,"Vascular Hydraulic Conductivity, Lp", Lp,"mm2 s/kg");
    
    // Interstitial permeability - mm3 s / kg
    kappa = 3.10e-05; //(1.7 / 0.1333) * 1e-5;
    outputf(ift,"Interstitial Permeability, kappa", kappa,"mm3 s/kg");
    
    // Far-field interstitial pressure - mmHg
    iP = 0;
    outputf(ift,"Far-field Interstitial Pressure", iP,"mmHg");
    
    // Oncotic pressure gradient
    PI_b = 20;
    PI_v = 15;
    sig = 0.82; // oncotic reflection coefficient
    oncotic = sig*(PI_b - PI_v);
    outputf(ift,"Vascular Oncotic Pressure", PI_b,"mmHg");
    outputf(ift,"Interstitial Oncotic Pressure", PI_v,"mmHg");
    outputf(ift,"Oncotic Reflection Coefficient", sig,"");
    
}
