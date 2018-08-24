//
//  netStatistics.cpp
//  Vascular Flows 2.0
//
//  Statistical analysis of inputted network file
//
//  Created by Paul Sweeney on 03/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>

#include "output_data.h"

extern string root;
extern float alx,aly,alz;
extern double pi,netVol,vasDen;
extern vec lseg,rseg;

void netStatistics(const string &name, const string &p1name, const int &p1Total, const vec &p1, const string &p2name, const int &p2Total, const vec &p2, const int &logScale, const string &p3name, const int &p3Total, const uvec &p3) {

    
    // Param. 1 statistics
    double p1max = max(p1);
    double p1min = min(p1);
    if (p1min == 0.0)    {
        cout<<"\t *** Warning: P1 minimum has been computed as 0um ***"<<endl;
    }
    double p1mean = mean(p1);
    double p1SD = stddev(p1);
    
    
    // Param. 2 statistics
    double p2max = max(p2);
    double p2min = min(p2);
    if (p2min == 0.0)    {
        cout<<"\t *** Warning: P2 minimum has been computed as 0um ***"<<endl;
    }
    double p2mean = mean(p2);
    double p2SD = stddev(p2);
    
    // Param. 3 statistics
    double p3max = max(p3);
    double p3min = min(p3);
    if (p3min == 0.0)    {
        cout<<"\t *** Warning: P3 minimum has been computed as 0um ***"<<endl;
    }
    double p3mean = mean(p3);
    vec temp2 = conv_to<vec>::from(p3);
    double p3SD = stddev(temp2);
    
    
    // Histogram Data
    int range = round(p1max)-round(p1min)+1;

    vec p1range = linspace<vec>(round(p1min),round(p1max),range);
    vec p1Histo = conv_to<vec>::from(hist(p1,p1range));

    range = round(p2max)-round(p2min)+1;
    vec p2range = linspace<vec>(round(p2min),round(p2max),range);
    vec p2Histo = conv_to<vec>::from(hist(p2,p2range));
    
    range = round(p3max)-round(p3min)+1;
    vec p3range = linspace<vec>(1,round(p3max),range);
    vec p3Histo = conv_to<vec>::from(p3);
    p3Histo = conv_to<vec>::from(hist(p3Histo,p3range));
    
    p1Histo /= p1Total;
    p2Histo /= p2Total;
    p3Histo /= p3Total;
    
    
    FILE *ofp1;
    
    string rootname = root + name;
    
    ofp1 = fopen(rootname.c_str(),"w");
    
    fprintf(ofp1,"Histogram Data\n");
    fprintf(ofp1,"%s data:\n",p1name.c_str());
    fprintf(ofp1,"Mean = %g , S.D. = %g , min = %g , max  = %g\n", p1mean,p1SD,p1min,p1max);
    fprintf(ofp1,"value  %% cumul. %%\n");
    for (int i = 0; i < p1range.n_elem; i++)  {
        fprintf(ofp1,"%f \t %f\n", p1range(i), p1Histo(i));
    }
    fprintf(ofp1,"\n%s data:\n",p2name.c_str());
    if (logScale == 1)  {
        fprintf(ofp1,"Mean = %g , S.D. = %g , min = %g , max  = %g\n", mean(exp(p2)),stddev(exp(p2)),min(exp(p2)),max(exp(p2)));
    }
    else{
        fprintf(ofp1,"Mean = %g , S.D. = %g , min = %g , max  = %g\n", p2mean,p2SD,p2min,p2max);
    }
    fprintf(ofp1,"value  %% cumul. %%\n");
    for (int i = 0; i < p2range.n_elem; i++) {
        fprintf(ofp1,"%f \t %f\n", p2range(i), p2Histo(i));
    }
    fprintf(ofp1,"\n%s data:\n",p3name.c_str());
    fprintf(ofp1,"Mean = %g , S.D. = %g , min = %g , max  = %g\n", p3mean,p3SD,p3min,p3max);
    fprintf(ofp1,"value  %% cumul. %%\n");
    for (int i = 0; i < p3range.n_elem; i++) {
        fprintf(ofp1,"%f \t %f\n", p3range(i), p3Histo(i));
    }
    
    fclose(ofp1);
    
}