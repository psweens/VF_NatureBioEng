//
//  bicgstab.cpp
//  Vascular-Flow
//
//  Biconjugate gradient stabilisation method taken from Greens O2 by Tim Secomb
//  and update for the armadillo library
//
//  Method has also been updated to be used for large symmetric matrices - l_bicgstab
//
//  Created by Paul Sweeney on 06/07/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <armadillo>
#include "global_variables.h"

using namespace arma;

double bicgstab(fmat &a, vec &b, vec &x, int &n, double &eps, int &itmax)     {
    
    double lu,lunew,beta,delta,er,gamma1,t1,t2,err=0.;
    int i,j,kk;
    vec r = zeros<vec>(n);
    vec rs = zeros<vec>(n);
    vec v = zeros<vec>(n);
    vec s = zeros<vec>(n);
    vec t = zeros<vec>(n);
    vec p = zeros<vec>(n);
    lu = 0.;
    for(i = 0; i < n; i++){
        r(i) = 0.;
        for(j = 0; j < n; j++)  {
            r(i) += a(i,j)*x(j);
        }
        r(i) -= b(i);
        p(i) = r(i);
        rs(i) = 1.;
        lu += r(i)*rs(i);
    }
    kk = 1;
    do  {
        t1 = 0.;
        for(i = 0; i < n; i++)  {
            v(i) = 0.;
            for(j = 0; j < n; j++)  {
                v(i) += a(i,j)*p(j);
            }
            t1 += v(i)*rs(i);
        }
        delta = -lu/t1;
        for(i = 0; i < n; i++) {
            s(i) = r(i) + delta*v(i);
        }
        for(i = 0; i < n; i++){
            t(i) = 0.;
            for(j = 0; j < n; j++)  {
                t(i) += a(i,j)*s(j);
            }
        }
        t1 = 0.;
        t2 = 0.;
        for(i = 0; i < n; i++){
            t1 += s(i)*t(i);
            t2 += t(i)*t(i);
        }
        gamma1 = -t1/t2;

        err = 0.;
        lunew = 0.;
        for(i = 0; i < n; i++){
            r(i) = s(i) + gamma1*t(i);
            er = delta*p(i) + gamma1*s(i);
            x(i) += er;
            if(abs(er) > err) {
                err = abs(er);
            }
            lunew += r(i)*rs(i);
        }
        beta = lunew*delta/(lu*gamma1);
        lu = lunew;
        
        for(i = 0; i < n; i++) {
            p(i) = r(i) + beta*(p(i)+gamma1*v(i));
        }
        printf("Loop %i: Error: %e\n",kk,err);
        kk += 1;
    }
    while(kk < itmax && err > eps);

    if(err > eps) printf("*** Warning: linear solution using BICGSTB not converged, err = %e\n",err);
    return err;
}


double bicgstab_sym(vec &b, vec &x, int &n, double &eps, int &itmax)   {
    
    double lu,lunew,beta,delta,er,gamma1,t1,t2,err=0.;
    int i,j,kk;
    vec r = zeros<vec>(n);
    vec rs = zeros<vec>(n);
    vec v = zeros<vec>(n);
    vec s = zeros<vec>(n);
    vec t = zeros<vec>(n);
    vec p = zeros<vec>(n);
    lu = 0.;
    cout<<"first loop"<<endl;
    for(i = 0; i < n; i++){
        
        FILE *ofp;
        char sim[64];
        sprintf(sim,"%i.txt",i);
        string rootname = "/Volumes/Samsung SSD/G_Data/";
        rootname += sim;
        ofp = fopen(rootname.c_str(),"r");
        
        double vector = 0.;
        double t_r = 0.;
        for (j = 0; j < n; j++) {
            fscanf(ofp,"%lf\n",&vector);
            t_r += vector*x(j);
        }
        fclose(ofp);
        t_r -= b(i);
        r(i) = t_r;
        p(i) = t_r;
        rs(i) = 1.;
        lu += t_r*rs(i);
        
    }
    kk = 1;
    do  {
        t1 = 0.;
        cout<<"second loop"<<endl;
        for(i = 0; i < n; i++)  {
            
            FILE *ofp;
            char sim[64];
            sprintf(sim,"%i.txt",i);
            string rootname = "/Volumes/Samsung SSD/G_Data/";
            rootname += sim;
            ofp = fopen(rootname.c_str(),"r");
            
            double vector = 0.;
            double t_v = 0.;
            for (j = 0; j < n; j++) {
                fscanf(ofp,"%lf\n",&vector);
                t_v += vector*p(j);
            }
            fclose(ofp);
            
            v(i) = t_v;
            t1 += t_v*rs(i);
        }
        delta = -lu/t1;
        for(i = 0; i < n; i++) {
            s(i) = r(i) + delta*v(i);
        }
        cout<<"third loop"<<endl;
        for(i = 0; i < n; i++){
            
            FILE *ofp;
            char sim[64];
            sprintf(sim,"%i.txt",i);
            string rootname = "/Volumes/Samsung SSD/G_Data/";
            rootname += sim;
            ofp = fopen(rootname.c_str(),"r");
            
            double vector = 0.;
            double t_t = 0.;
            for (j = 0; j < n; j++) {
                fscanf(ofp,"%lf\n",&vector);
                t_t += vector*s(j);
            }
            fclose(ofp);
            
            t(i) = t_t;
        }
        t1 = 0.;
        t2 = 0.;
        for(i = 0; i < n; i++){
            t1 += s(i)*t(i);
            t2 += t(i)*t(i);
        }
        gamma1 = -t1/t2;
        err = 0.;
        lunew = 0.;
        for(i = 0; i < n; i++){
            r(i) = s(i) + gamma1*t(i);
            er = delta*p(i) + gamma1*s(i);
            x(i) += er;
            if(fabs(er) > err) err = fabs(er);
            lunew += r(i)*rs(i);
        }
        beta = lunew*delta/(lu*gamma1);
        lu = lunew;
        for(i = 0; i < n; i++) {
            p(i) = r(i) + beta*(p(i)+gamma1*v(i));
        }
        printf("Loop %i: Error: %e\n",kk,err);
        kk += 1;
    }
    while(kk < itmax && err > eps);
    
    if(err > eps) printf("*** Warning: linear solution using BICGSTB not converged, err = %e\n",err);
    return err;
}


double l_bicgstab(vec &b, vec &x, int &n, double &eps, int &itmax)   {
    
    double lu,lunew,beta,delta,er,gamma1,t1,t2,err=0.;
    int i,j,kk;
    vec r = zeros<vec>(n);
    vec rs = zeros<vec>(n);
    vec v = zeros<vec>(n);
    vec s = zeros<vec>(n);
    vec t = zeros<vec>(n);
    vec p = zeros<vec>(n);
    lu = 0.;
    cout<<"first loop"<<endl;
    for(i = 0; i < n; i++){
        
        FILE *ofp;
        char sim[64];
        sprintf(sim,"%i.txt",i);
        string rootname = "/Volumes/Samsung SSD/G_Data/";
        rootname += sim;
        ofp = fopen(rootname.c_str(),"r");
        
        double vector = 0.;
        double t_r = 0.;
        for (j = 0; j < n; j++) {
            fscanf(ofp,"%lf\n",&vector);
            t_r += vector*x(j);
        }
        fclose(ofp);
        t_r -= b(i);
        r(i) = t_r;
        p(i) = t_r;
        rs(i) = 1.;
        lu += t_r*rs(i);
        
    }
    kk = 1;
    do  {
        t1 = 0.;
        cout<<"second loop"<<endl;
        for(i = 0; i < n; i++)  {
            
            FILE *ofp;
            char sim[64];
            sprintf(sim,"%i.txt",i);
            string rootname = "/Volumes/Samsung SSD/G_Data/";
            rootname += sim;
            ofp = fopen(rootname.c_str(),"r");
            
            double vector = 0.;
            double t_v = 0.;
            for (j = 0; j < n; j++) {
                fscanf(ofp,"%lf\n",&vector);
                t_v += vector*p(j);
            }
            fclose(ofp);
            
            v(i) = t_v;
            t1 += t_v*rs(i);
        }
        delta = -lu/t1;
        for(i = 0; i < n; i++) {
            s(i) = r(i) + delta*v(i);
        }
        cout<<"third loop"<<endl;
        for(i = 0; i < n; i++){
            
            FILE *ofp;
            char sim[64];
            sprintf(sim,"%i.txt",i);
            string rootname = "/Volumes/Samsung SSD/G_Data/";
            rootname += sim;
            ofp = fopen(rootname.c_str(),"r");
            
            double vector = 0.;
            double t_t = 0.;
            for (j = 0; j < n; j++) {
                fscanf(ofp,"%lf\n",&vector);
                t_t += vector*s(j);
            }
            fclose(ofp);
            
            t(i) = t_t;
        }
        t1 = 0.;
        t2 = 0.;
        for(i = 0; i < n; i++){
            t1 += s(i)*t(i);
            t2 += t(i)*t(i);
        }
        gamma1 = -t1/t2;
        err = 0.;
        lunew = 0.;
        for(i = 0; i < n; i++){
            r(i) = s(i) + gamma1*t(i);
            er = delta*p(i) + gamma1*s(i);
            x(i) += er;
            if(fabs(er) > err) err = fabs(er);
            lunew += r(i)*rs(i);
        }
        beta = lunew*delta/(lu*gamma1);
        lu = lunew;
        for(i = 0; i < n; i++) {
            p(i) = r(i) + beta*(p(i)+gamma1*v(i));
        }
        printf("Loop %i: Error: %e\n",kk,err);
        kk += 1;
    }
    while(kk < itmax && err > eps);
    
    if(err > eps) printf("*** Warning: linear solution using BICGSTB not converged, err = %e\n",err);
    return err;
    
}
