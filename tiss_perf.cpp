//
//  tiss_perf.cpp
//  Vascular-Flow
//
//  Created by Paul Sweeney on 24/04/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "inter_var.h"
#include "global_variables.h"
#include "global_params.h"
#include "output_data.h"
#include "inter_func.h"
#include "misc_func.h"

void tiss_perf()    {
    
    double rho = 1.05*1e-12;
    double discrete_inflow = accu(BCflow(find(BCflow > 0.0)));
    
    if (nseg < 1e4) {
        
        double val = 20.; // length of discretised squares
        double disc_x = floor(alx / val);
        double disc_y = floor(aly / val);
        double disc_z = floor(alz / val);
        double xl = alx / disc_x;
        double yl = aly / disc_y;
        double zl = alz / disc_z;
        
        int pnts = 2*(disc_x*disc_y + disc_x*disc_z + disc_y*disc_z);
        mat pnode = zeros<mat>(3,pnts);
        int cntr = 0;
        // x,y faces
        ivec face = zeros<ivec>(pnts);
        for (int i = 0; i < disc_x; i++)    {
            for (int j = 0; j < disc_y; j++)    {
                pnode(0,cntr) = (i + 0.5)*xl;
                pnode(1,cntr) = (j + 0.5)*yl;
                face(cntr) = 1;
                cntr += 1;
                pnode(0,cntr) = (i + 0.5)*xl;
                pnode(1,cntr) = (j + 0.5)*yl;
                pnode(2,cntr) = alz;
                face(cntr) = 1;
                cntr += 1;
            }
        }
        // y,z faces
        for (int i = 0; i < disc_z; i++)    {
            for (int j = 0; j < disc_y; j++)    {
                pnode(1,cntr) = (j + 0.5)*yl;
                pnode(2,cntr) = (i + 0.5)*zl;
                face(cntr) = 2;
                cntr += 1;
                pnode(0,cntr) = alx;
                pnode(1,cntr) = (j + 0.5)*yl;
                pnode(2,cntr) = (i + 0.5)*zl;
                face(cntr) = 2;
                cntr += 1;
            }
        }
        // x,z faces
        for (int i = 0; i < disc_x; i++)    {
            for (int j = 0; j < disc_z; j++)    {
                pnode(0,cntr) = (i + 0.5)*xl;
                pnode(2,cntr) = (j + 0.5)*zl;
                face(cntr) = 3;
                cntr += 1;
                pnode(0,cntr) = (i + 0.5)*xl;
                pnode(1,cntr) = aly;
                pnode(2,cntr) = (j + 0.5)*zl;
                face(cntr) = 3;
                cntr += 1;
            }
        }
        
        vec surfv = zeros<vec>(pnts);
        for (int i = 0; i < pnts; i++)  {
            vec x = pnode.col(i);
            surfv(i) = DCeval(x,1);
        }
        
        // Compute average flow at each face
        double interstitial_inflow = 0.0;
        for (int i = 0; i < pnts; i++)  {
            if (surfv(i) > 0.0) {
                if (face(i) == 1)    {
                    interstitial_inflow += (xl * yl * surfv(i));
                }
                else if (face(i) == 2)    {
                    interstitial_inflow += (yl * zl * surfv(i));
                }
                else if (face(i) == 3)    {
                    interstitial_inflow += (xl * zl * surfv(i));
                }
            }
        }
        interstitial_inflow *= (60/1e6); // um3/s --> nl/min
        double m = alx*aly*alz*rho;
        double tissue_perf = (discrete_inflow + interstitial_inflow)*(100/m)*1e-6;
        
        net2amira("Interstitial/Nodes2Amira.txt","VessTyp", pnts, pnts, pnode, zeros<uvec>(pnts), zeros<uvec>(pnts), zeros<vec>(pnts), zeros<vec>(pnts));
        
    }
    else {
        
        // Inflow from boundary nodes
        ivec flag = zeros<ivec>(nnodbc);
        mat perf_xyz;
        string name = "Read_Hull";
        hull_read(name,flag,perf_xyz,0);
        inflow = accu(BCflow(find(BCflow > 0. && flag == 1)));
        outputf(ift,"Vascular Inflow", inflow,"nl/min");
        
        // Inflow along surface of tumour
        name = "Read_Hull";
        hull_read(name,flag,perf_xyz,1);
        netVol = sum(lseg % (diam % diam))*M_PI/4;
        vasDen = 100*netVol/vol;
        outputf(ift,"Tissue Volume", vol*1e-9,"mm3");
        outputf(ift,"(Updated) Vascular Density", vasDen,"\%");
        int elem = (int) perf_xyz.n_cols / 2;
        mat hull = perf_xyz.cols(0,elem-1);
        perf_xyz = perf_xyz.cols(elem,perf_xyz.n_cols-1);
        vec radii = zeros<vec>(perf_xyz.n_cols);
        flag = zeros<ivec>(perf_xyz.n_cols);
        radii.fill(1);
        ivec flag_circle = sphere_packing(1.05, radii, perf_xyz);
        outputf(ift,"Min. Boundary Sphere", min(radii),"um");
        outputf(ift,"Max. Boundary Sphere", max(radii),"um");
        
        vec surfv = zeros<vec>(perf_xyz.n_cols);
        
        for (int inod = 0; inod < perf_xyz.n_cols; inod++)  {
            vec x = perf_xyz.col(inod);
            vec xx = hull.col(inod);
            double press1 = DCeval(x,0);
            double press2 = DCeval(xx,0);
            if (press1 > press2)    {
                surfv(inod) = DCeval(x,1);
                flag(inod) = 1;
            }
        }
        
        double m = vol*rho; // Calculate mass
        double int_inflow = accu(surfv(find(flag == 1)) % (M_PI*pow(radii(find(flag == 1)),2)))*(60/1e6); // Assumed contributing cross-sectional area
        outputf(ift,"Interstitial Inflow", int_inflow,"nl/min");
        inflow += int_inflow;
        RBF = inflow*(100/m)*1e-6;
        
        outputf(ift,"Tissue Perfusion", RBF,"ml/min/100g");
        
        
        FILE *ipf;
        
        string filename = "Interstitial/inflow_xyz_radii.txt";
        filename = root + filename;
        
        ipf = fopen(filename.c_str(),"w");
        
        for (int is = 0; is < perf_xyz.n_cols; is++) {
            fprintf(ipf,"%f %f %f %f\n",perf_xyz(0,is),perf_xyz(1,is),perf_xyz(2,is),radii(is));
        }
        
        fclose(ipf);
    }
    
}
