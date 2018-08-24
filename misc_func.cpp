//
//  pearson_coeff.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 31/08/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <sstream>

#include "global_variables.h"


using namespace std;
using namespace arma;


int detect_col(FILE *ifp)    {
    
    int max=200;
    char bb[200];
    string s;
    fpos_t position;
    fgetpos(ifp,&position);
    fgets(bb,max,ifp);
    s = bb;
    istringstream is( s );
    double n;
    int cnt = 0;
    while( is >> n ) {
        cnt += 1;
    }
    fsetpos(ifp,&position);
    
    return cnt;
    
}


void outputf(FILE *ofp, const string &text, const double &num, const string &unit)   {
    
    string rootname = root + "Interstitial/GreensLog.txt";
    ofp = fopen(rootname.c_str(),"a");
    
    if (round(num) == num)  {
        printf("%s : %i %s\n",text.c_str(),(int) num,unit.c_str());
        fprintf(ofp,"%s : %i %s\n",text.c_str(),(int) num,unit.c_str());
    }
    else if (num < 1e-3)    {
        printf("%s : %.2e %s\n",text.c_str(),num,unit.c_str());
        fprintf(ofp,"%s : %.2e %s\n",text.c_str(),num,unit.c_str());
    }
    else {
        printf("%s : %.2f %s\n",text.c_str(),num,unit.c_str());
        fprintf(ofp,"%s : %.2f %s\n",text.c_str(),num,unit.c_str());
    }
    
    fclose(ofp);

}

void init_log(FILE *ofp)    {
    
    string rootname = root + "Interstitial/GreensLog.txt";
    ofp = fopen(rootname.c_str(),"w");
    fclose(ofp);
    
}

int max_double(const double &a, const double &b) {
    
    if (a > b)  {
        return 1;
    }
    else if (a < b) {
        return 2;
    }
    else {
        return 1;
    }
    
}


ivec sphere_packing(const double grow_decay, vec &r0, const mat &snode)   {
    
    int circ_count = 0;
    int ns = (int) snode.n_cols;
    ivec flag_sphere = zeros<ivec>(ns);
    vec old_r0 = zeros<vec>(ns);
    flag_sphere = zeros<ivec>(ns);
    
    int old_circ_count = -1;
    int cnt = 0;
    while (circ_count < ns)  {
        
        double dist = 0.0;
        for (int is = 0; is < ns; is++) {
            if (flag_sphere(is) == 0)   {
                int cntr = 0;
                for (int js = 0; js < ns; js++) {
                    if (is != js && flag_sphere(js) == 0 && grow_decay < 1.)   {
                        dist = sqrt(pow(snode(0,is)-snode(0,js),2) + pow(snode(1,is)-snode(1,js),2) + pow(snode(2,is)-snode(2,js),2)) - r0(is) - r0(js);
                        if (dist < 0.0)  {
                            cntr += 1;
                        }
                    }
                    else if (is != js && grow_decay > 1.)   {
                        dist = sqrt(pow(snode(0,is)-snode(0,js),2) + pow(snode(1,is)-snode(1,js),2) + pow(snode(2,is)-snode(2,js),2)) - r0(is) - r0(js);
                        if (dist < 0.0)  {
                            cntr += 1;
                        }
                    }
                }
                if (cntr == 0 && grow_decay < 1.)  {
                    flag_sphere(is) = 1;
                }
                else if (cntr == 0 && grow_decay > 1.)  {
                    flag_sphere(is) = 0;
                }
                else if (cntr > 0 && grow_decay > 1.)   {
                    flag_sphere(is) = 1;
                    r0(is) *= (2 - grow_decay);
                }
                if (cnt > 10 && grow_decay < 1.)   {
                    flag_sphere(is) = 2;
                }
                else if (cnt > 1e2 && grow_decay > 1.)  {
                    flag_sphere(is) = 2;
                }
            }
        }
        
        old_r0 = r0;
        r0(find(flag_sphere == 0)) *= grow_decay;
        circ_count = (int) accu(flag_sphere);
        
        if (old_circ_count == circ_count)   {
            cnt += 1;
        }
        else {
            cnt = 0;
        }
        old_circ_count = circ_count;
    }
    
    
    return flag_sphere;
}

void time_check(FILE *ift, double &run_start, const string &text)   {
    
    double time_span = (clock()-run_start)/(double) CLOCKS_PER_SEC;
    if (time_span > 3600)    {
        outputf(ift,text, time_span/3600,"hours");
    }
    else if (time_span > 60)    {
        outputf(ift,text, time_span/60,"minutes");
    }
    else    {
        outputf(ift,text, time_span,"seconds");
    }
    
}
