//
//  finite_mesh.cpp
//  Vascular-Flow
//
//  Code to calculate perfusion across a mesh to be compared with MRI data
//
//  Created by Paul Sweeney on 22/06/2017.
//  Copyright © 2017 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "inter_var.h"
#include "inter_func.h"
#include "global_params.h"
#include "global_variables.h"
#include "misc_func.h"
#include "discrete_flow.h"

void finite_mesh(const int &disc_x, const int &disc_y, const int &disc_z)  {
    
    double x_length = alx/disc_x;
    double y_length = aly/disc_y;
    double z_length = alz/disc_z;
    
    cube pressure = zeros<cube>(disc_x+1,disc_y+1,disc_z+1);
    cube velocity = zeros<cube>(disc_x+1,disc_y+1,disc_z+1);
    cube perfusion = zeros<cube>(disc_x,disc_y,disc_z);
    
    double xmax = 0.;
    for (int x = 0; x <= disc_x; x++)    {
        
        double ymax = 0.;
        for (int y = 0; y <= disc_y; y++)    {
            
            double zmax = 0.;
            for (int z = 0; z <= disc_z; z++)   {
                
                vec coord = zeros<vec>(3);
                coord(0) = xmax;
                coord(1) = ymax;
                coord(2) = zmax;
                pressure(x,y,z) = DCeval(coord,0) + iP;
                velocity(x,y,z) = DCeval(coord,1)*1e3;
                
                zmax += z_length;
            }
            
            ymax += y_length;
        }
        
        xmax += x_length;
    }

    // Cube surface areas
    double xy_surf = x_length * y_length;
    double xz_surf = x_length * z_length;
    double yz_surf = y_length * z_length;
    
    cube streamline_u = zeros<cube>(disc_x,disc_y,disc_z);
    cube streamline_v = zeros<cube>(disc_x,disc_y,disc_z);
    cube streamline_w = zeros<cube>(disc_x,disc_y,disc_z);
    for (int x = 0; x < disc_x; x++)   {
        for (int y = 0; y < disc_y; y++)    {
            for (int z = 0; z < disc_z; z++)    {
                
                int flag = max_double(pressure(x,y,z),pressure(x+1,y,z));
                if (flag == 1) {
                    perfusion(x,y,z) += (yz_surf * velocity(x,y,z));
                }
                else {
                    perfusion(x,y,z) += (yz_surf * velocity(x+1,y,z));
                }
                
                flag = max_double(pressure(x,y,z),pressure(x,y+1,z));
                if (flag == 1) {
                    perfusion(x,y,z) += (xz_surf * velocity(x,y,z));
                }
                else {
                    perfusion(x,y,z) += (xz_surf * velocity(x,y+1,z));
                }
                
                flag = max_double(pressure(x,y,z),pressure(x,y,z+1));
                if (flag == 1) {
                    perfusion(x,y,z) += (xy_surf * velocity(x,y,z));
                }
                else {
                    perfusion(x,y,z) += (xy_surf * velocity(x,y,z+1));
                }
                
                streamline_u(x,y,z) = pressure(x,y,z) - pressure(x+1,y,z);
                streamline_v(x,y,z) = pressure(x,y,z) - pressure(x,y+1,z);
                streamline_w(x,y,z) = pressure(x,y,z) - pressure(x,y,z+1);
                
            }
        }
    }
    
    double vol = xy_surf * xz_surf * yz_surf; // Cube volume - um3
    double rho = 1.05*1e-12; // Tissue density - g / um3
    double mass = vol*rho;
    
    perfusion = perfusion*(100*60/mass)*1e-12;
    outputf(ift,"Min. Interstitial Perfusion", perfusion.min(),"ml/min/100g");
    outputf(ift,"Max. Interstitial Perfusion", perfusion.max(),"ml/min/100g");
    
    
    // Include vascular contribution to perfusion
    //cube vasc_perf = discretePerfusion(disc_x, disc_y, disc_z);
    //perfusion = perfusion + vasc_perf;
    
    
    string slice_folder = root + "Interstitial/Slices/Perfusion/";
    string slice_u_folder = root + "Interstitial/Slices/Streamline/u/";
    string slice_v_folder = root + "Interstitial/Slices/Streamline/v/";
    string slice_w_folder = root + "Interstitial/Slices/Streamline/w/";
    
    char slice[64];
    for (int z = 0; z < disc_z; z++)    {
        
        sprintf(slice,"%i.txt",z);
        string slice_name = slice_folder + slice;
        FILE *ofp = fopen(slice_name.c_str(),"w");
        for (int x = 0; x < disc_x; x++)  {
            for (int y = 0; y < disc_y; y++)    {
                fprintf(ofp,"%e\t",perfusion(x,y,z));
            }
            fprintf(ofp,"\n");
        }
        fclose(ofp);
        
        
        sprintf(slice,"%i.txt",z);
        string slice_u = slice_u_folder + slice;
        FILE *osu = fopen(slice_u.c_str(),"w");
        for (int x = 0; x < disc_x; x++)  {
            for (int y = 0; y < disc_y; y++)    {
                fprintf(osu,"%e\t",streamline_u(x,y,z));
            }
            fprintf(osu,"\n");
        }
        fclose(osu);
        
        
        sprintf(slice,"%i.txt",z);
        string slice_v = slice_v_folder + slice;
        FILE *osv = fopen(slice_v.c_str(),"w");
        for (int x = 0; x < disc_x; x++)  {
            for (int y = 0; y < disc_y; y++)    {
                fprintf(osv,"%e\t",streamline_v(x,y,z));
            }
            fprintf(osv,"\n");
        }
        fclose(osv);
        
        
        sprintf(slice,"%i.txt",z);
        string slice_w = slice_w_folder + slice;
        FILE *osw = fopen(slice_w.c_str(),"w");
        for (int x = 0; x < disc_x; x++)  {
            for (int y = 0; y < disc_y; y++)    {
                fprintf(osw,"%e\t",streamline_w(x,y,z));
            }
            fprintf(osw,"\n");
        }
        fclose(osw);
        
    }
    
}
