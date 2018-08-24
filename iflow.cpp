//
//  iflow.cpp
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
#include "output_data.h"
#include "misc_func.h"

void iflow()    {
    
    clock_t sub_run_start;
    sub_run_start = clock();

    float alpha = 0.1333;
    
    s_p = zeros<vec>(ns);
    vec Jv = zeros<vec>(ns);
    s_q = zeros<vec>(ns);

    double qerror = 1e3;
    double qerror_old;
    int cnt_err = 5;
    int cntr = 1;
    
    FILE *ofp;
    string filename = "Interstitial/Convergence_Data/Convergence_Data_0.txt";
    filename = root + filename;
    ofp = fopen(filename.c_str(),"w");
    for (int is = 0; is < ns; is++) {
        fprintf(ofp,"%f %f %e %f\n",s_segpress(is), Jv(is), s_q(is), old_q(is));
    }
    fclose(ofp);
    

    
    if (ns < 5e4)   {
        G = zeros<mat>(ns,ns);
        double t_r0 = 0.;
        vec t_snode = zeros<vec>(3);
        for (int i = 0; i < ns; i++)    {
            t_r0 = r0(i);
            t_snode = snode.col(i);
            for (int j = 0; j < ns; j++)    {
                G(i,j) = eval_G(t_r0, kappa, t_snode, snode.col(j)) - (2*M_PI * (kappa / Lp) * eval_deltaG(t_r0, kappa, t_snode, snode.col(j))) * 1e-6;
            }
        }
    }
    else if (ns < 1e5)  {
        float_G = zeros<fmat>(ns,ns);
        double t_r0 = 0.;
        vec t_snode = zeros<vec>(3);
        for (int i = 0; i < ns; i++)    {
            t_r0 = r0(i);
            t_snode = snode.col(i);
            for (int j = 0; j < ns; j++)    {
                float_G(i,j) = eval_G(t_r0, kappa, t_snode, snode.col(j)) - (2*M_PI * (kappa / Lp) * eval_deltaG(t_r0, kappa, t_snode, snode.col(j))) * 1e-6;
            }
        }
    }
    else {
        cout<<"start G"<<endl;
        double t_r0 = 0.;
        vec t_snode = zeros<vec>(3);
        for (int i = 0; i < ns; i++)    {
            FILE *ofp;
            char sim[64];
            sprintf(sim,"%i.txt",i);
            string rootname = "/Volumes/Samsung SSD/G_Data/";
            rootname += sim;
            ofp = fopen(rootname.c_str(),"w");
            
            t_r0 = r0(i);
            t_snode = snode.col(i);
            
            for (int j = 0; j < ns; j++)    {
                
                double temp_G = eval_G(t_r0, kappa, t_snode, snode.col(j)) - (2*M_PI * (kappa / Lp) * eval_deltaG(t_r0, kappa, t_snode, snode.col(j)) * 1e-6);
                
                fprintf(ofp,"%f\n",temp_G);
                
            }
            fclose(ofp);
        }
        
        
    }
    vec rhs = (s_segpress - iP - oncotic)*alpha;
    
    while (abs(qerror) > 1e-8 && cnt_err <= 5) {
        
        FILE *ifp;
        filename = "Interstitial/Convergence_Data/Convergence_Data_";
        char sim[2];
        sprintf(sim,"%i",cntr);
        filename += sim;
        filename += ".txt";
        filename = root + filename;
        
        ifp = fopen(filename.c_str(),"w");
        
        
        
        
        s_q.zeros();
        if (ns < 5e4)   {
            s_q = solve(G,rhs);
            s_p = iP + G*s_q;
        }
        else if (ns < 1e5)  {
            int itmax = 5;
            double eps = 1e-8;
            bicgstab(float_G, rhs, s_q, ns, eps, itmax);
            s_p = iP + float_G*s_q;
        }
        else {

            int itmax = 3;
            double eps = 1e-8;
            cout<<"start solver"<<endl;
            l_bicgstab(rhs, s_q, ns, eps, itmax);
            
            cout<<"system solved"<<endl;
            for (int i = 0; i < ns; i++)  {
                FILE *ofp;
                char sim[64];
                sprintf(sim,"%i.txt",i);
                string rootname = "/Volumes/Samsung SSD/G_Data/";
                rootname += sim;
                ofp = fopen(rootname.c_str(),"r");
                double temp_G = 0.;
                
                for (int j = i; j < ns; j++)  {
                    fscanf(ofp,"%lf\n",&temp_G);
                    if (i == j) {
                        s_p(i) += temp_G*s_q(j);
                    }
                    else {
                        s_p(i) += 2*temp_G*s_q(j);
                    }
                }
                fclose(ofp);
                s_p(i) += iP;
            }
            
        }
        
        
        for (int is = 0; is < ns; is++) {
            if (s_p(is) < 0.0)  {
                printf("\n\t\t*** Warning : Negative Pressure Detected ***\n");
                is = ns;
            }
        }
        
        
        vec s_flow = zeros<vec>(ns);
        vec s_press = zeros<vec>(ns);
        for (int is = 0; is < ns; is++) {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (sname(is) == segname(iseg)) {
                    s_flow(is) = qq(iseg); // nl/min to um3/s
                    s_press(is) = segpress(iseg);
                }
            }
        }

        
        for (int i = 0; i < ns; i++)    {
            double t_r0 = r0(i);
            double t_rseg = s_rseg(i)*1e-3;
            double t_lseg = s_lseg(i)*1e-3;
            vec t_snode = snode.col(i);
            double t_Jv = 0.;
            for (int j = 0; j < ns; j++)    {
                t_Jv -= 2*M_PI*kappa*eval_deltaG(t_r0, kappa, t_snode, snode.col(j))* t_rseg * t_lseg * s_q(j)*1e-6;
            }
            Jv(i) = t_Jv;
        }
        
        
        s_press = s_press - (iK%Jv)/alpha;

        
        qerror_old = qerror;
        qerror = max(old_q - (old_q - Jv/Gamma));
        old_q -= Jv/Gamma;
        
        if (qerror == qerror_old)   {
            //cnt_err += 1;
        }
        cnt_err += 1;
        
        cout<<"Max Flow error = "<<qerror<<endl;
        
        for (int is = 0; is < ns; is++) {
            fprintf(ifp,"%f %e %e %e\n",s_press(is), Jv(is), s_q(is)/Gamma, old_q(is));
        }
        
        fclose(ifp);
        cntr += 1;
        
        rhs = (s_press - iP - oncotic)*alpha;
        
        // Output source coordinates, pressures and strengths
        FILE *ipf;
        
        filename = "Interstitial/Source_xyz_Press_Strength";
        filename += sim;
        filename += ".txt";
        filename = root + filename;
        
        ipf = fopen(filename.c_str(),"w");
        for (int is = 0; is < ns; is++) {
            fprintf(ipf,"%lli %f %f %f %f %f %f %e\n",sname(is),snode(0,is),snode(1,is),snode(2,is),r0(is),s_lseg(is),s_p(is),s_q(is));
        }
        
        fclose(ipf);
        
    }

    
}
