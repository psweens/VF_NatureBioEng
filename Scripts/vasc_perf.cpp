//
//  discretePerfusion.cpp
//  Vascular Flows
//
//  vasc_perf - calculates vascular contribution to tissue perfusion 
//
//  discretePerfusion - creates coarse MRI-style perfusion data (similar to 'finite_mesh')
//  Discretises tissue domain into a user defined mesh and calculates vascular perfusion for each cell
//  Matrices can be loaded into matlab to create videos. Needs to be adjusted for use in Amira.
//
//  Created by Paul Sweeney on 24/03/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>

#include "global_variables.h"
#include "global_params.h"
#include "discrete_flow.h"
#include "output_data.h"

void vasc_perf(const int &tag)    {
    
    
    // Calculate boundary flow
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (bcnod(inodbc) == ista(iseg)) {
                if (nodpress(ista(iseg)) > nodpress(iend(iseg)))  {
                    BCflow(inodbc) = qq(iseg);
                }
                else if (nodpress(ista(iseg)) < nodpress(iend(iseg))) {
                    BCflow(inodbc) = -qq(iseg);
                }
            }
            else if (bcnod(inodbc) == iend(iseg))   {
                if (nodpress(ista(iseg)) < nodpress(iend(iseg)))  {
                    BCflow(inodbc) = qq(iseg);
                }
                else if (nodpress(ista(iseg)) > nodpress(iend(iseg))) {
                    BCflow(inodbc) = -qq(iseg);
                }
            }
        }
    }
    
    
    rho = 1.05*1e-12; // Density of tissue - g/ml
    double m = 0.; // Mass of tissue - defined below
    
    // 0 - Calculates vascular perfusion using inflow from all boundary nodes
    // 1 - Calculates vascular perfusion using boundary nodes located in a square surface of user-defined thickness around the network (for cuboid networks)
    // 2 - Calculates vascular perfusion using boundary nodes located on a pre-defined, network-specific convex hull of the network
    if (tag == 0)   {
        
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (bcnod(inodbc) == ista(iseg)) {
                    if (nodpress(ista(iseg)) > nodpress(iend(iseg)))  {
                        BCflow(inodbc) = qq(iseg);
                    }
                    else if (nodpress(ista(iseg)) < nodpress(iend(iseg))) {
                        BCflow(inodbc) = -qq(iseg);
                    }
                }
                else if (bcnod(inodbc) == iend(iseg))   {
                    if (nodpress(ista(iseg)) < nodpress(iend(iseg)))  {
                        BCflow(inodbc) = qq(iseg);
                    }
                    else if (nodpress(ista(iseg)) > nodpress(iend(iseg))) {
                        BCflow(inodbc) = -qq(iseg);
                    }
                }
            }
        }
        
        inflow = accu(BCflow(find(BCflow > 0.)));
        m = alx*aly*alz*rho;
        
    }
    else if (tag == 1)  {
        
        RBF = 0.;
        ivec flag = zeros<ivec>(nnodbc);
        ivec nod = zeros<ivec>(nnod);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (cnode(2,ista(iseg)) < 375. || cnode(2,iend(iseg)) < 375.)   {
                nod(ista(iseg)) += 1;
                nod(iend(iseg)) += 1;
            }
        }
        double margin = 0.1; // i.e. 10% margin rel. to tissue block
        double x1 = (1-margin)*alx;
        double x2 = margin*alx;
        double y1 = (1-margin)*aly;
        double y2 = margin*aly;
        double z1 = (1-margin)*alz;
        double z2 = margin*alz;
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            if ((cnode(0,bcnod(inodbc)) > x1 || cnode(0,bcnod(inodbc)) < x2 || cnode(1,bcnod(inodbc)) > y1 || cnode(1,bcnod(inodbc)) < y2 || cnode(2,bcnod(inodbc)) < z2) && cnode(2,bcnod(inodbc)) < z1) {
                flag(inodbc) = 1;
            }
        }
        vec flow = zeros<vec>(nseg);
        for (int iseg = 0; iseg < nseg; iseg++) {
            for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                if (flag(inodbc) == 1 && (ista(iseg) == bcnod(inodbc) || iend(iseg) == bcnod(inodbc)))   {
                    flow(iseg) = BCflow(inodbc);
                }
            }
        }
        for (int inod = 0; inod < nnod; inod++) {
            if (nod(inod) == 1 && nodtyp(inod) != 1)    {
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (inod == ista(iseg) && cnode(2,iend(iseg)) < z1) {
                        if (nodpress(ista(iseg)) > nodpress(iend(iseg)))  {
                            flow(inod) = qq(iseg);
                        }
                    }
                    else if (inod == iend(iseg) && cnode(2,ista(iseg)) < z1)   {
                        if (nodpress(ista(iseg)) < nodpress(iend(iseg)))  {
                            flow(inod) = qq(iseg);
                        }
                    }
                }
            }
        }
        
        vec sflag = zeros<vec>(nseg);
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            for (int iseg = 0; iseg < nseg; iseg++) {
                if (flag(inodbc) == 1 && (ista(iseg) == bcnod(inodbc) || iend(iseg) == bcnod(inodbc))) {
                    sflag(iseg) = 1;
                }
                else if (nod(ista(iseg)) == 1 && nodtyp(ista(iseg)) != 1)  {
                    sflag(iseg) = 1;
                }
                else if (nod(iend(iseg)) == 1 && nodtyp(iend(iseg)) != 1)  {
                    sflag(iseg) = 1;
                }
            }
        }
        
        inflow = accu(flow);
        m = alx*aly*alz*rho;
        
        double internal_flow = accu(BCflow(find(flag == 0)));
        cout<<"\t\t\t\tInternal Boundary Flow Contribution: "<<internal_flow<<" nl/min"<<endl;
        
    }
    else if (tag == 2)  {
        
        ivec flag = zeros<ivec>(nnodbc);
        mat temp;
        hull_read("Read_Hull",flag,temp,0);
        hull_read("Read_Hull",flag,temp,1);
        m = vol*rho;
        inflow = accu(BCflow(find(BCflow > 0. && flag == 1)));
        
        double internal_flow = accu(BCflow(find(flag == 0)));
        cout<<"\t\t\t\tInternal Boundary Flow Contribution: "<<internal_flow<<" nl/min"<<endl;
        
    }
    
    RBF = inflow*(100/m)*1e-6;
 
    cout<<"\t\t\t\tVascular Perfusion: "<<RBF<<" ml/min/100g"<<endl;
    
}


cube discretePerfusion(const int &disc_x, const int &disc_y, const int &disc_z)    {
    
    double x_length = alx/disc_x;
    double y_length = aly/disc_y;
    double z_length = alz/disc_z;
    
    double xmax = 0.;
    double ymax = 0.;
    double zmax = 0.;
    
    double rho = 1.05*1e-12;
    
    cube block_dens = zeros<cube>(disc_x, disc_y, disc_z);
    cube BC_dens = zeros<cube>(disc_x, disc_y, disc_z);
    
    cube disc_perf = zeros<cube>(disc_x, disc_y, disc_z);
    
    cube x_field = zeros<cube>(disc_x, disc_y, disc_z);
    cube y_field = zeros<cube>(disc_x, disc_y, disc_z);
    cube z_field = zeros<cube>(disc_x, disc_y, disc_z);
    
    double m = 0.0;
    
    double xmin, ymin, zmin;
    for (int x = 0; x < disc_x; x++)    {
        xmax += x_length;
        xmin = (xmax - x_length);
        
        ymax = 0.;
        for (int y = 0; y < disc_y; y++)    {
            ymax += y_length;
            ymin = (ymax - y_length);
            
            zmax = 0.;
            for (int z = 0; z < disc_z; z++)   {
                zmax += z_length;
                zmin = (zmax - z_length);
                
                
                
                ivec flagNod = zeros<ivec>(nnod);
                for (int inod = 0; inod < nnod; inod ++)    {
                    if (cnode(0,inod) >= xmin && cnode(0,inod) <= xmax && cnode(1,inod) >= ymin && cnode(1,inod) <= ymax && cnode(2,inod) >= zmin && cnode(2,inod) <= zmax)    {
                        flagNod(inod) = 1;
                    }
                }
                
                ivec flagSeg = zeros<ivec>(nseg);
                for (int iseg = 0; iseg < nseg; iseg ++)    {
                    if (flagNod(ista(iseg)) == 1 && flagNod(iend(iseg)) == 1)    {
                        flagSeg(iseg) = 1;
                    }
                }
                
                ivec sub_nodtyp = zeros<ivec>(nnod);
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (flagSeg(iseg) == 1) {
                        int inod1 = (int) ista(iseg);
                        int inod2 = (int) iend(iseg);
                        sub_nodtyp(inod1) += 1;
                        sub_nodtyp(inod2) += 1;
                    }
                }
                
                int cntr = 0;
                for (int iseg = 0; iseg < nseg; iseg++) {
                    if (sub_nodtyp(ista(iseg)) == 1)    {
                        if (nodpress(ista(iseg)) > nodpress(iend(iseg)))    {
                            disc_perf(x,y,z) += qq(iseg);
                            cntr += 1;
                        }
                    }
                    else if (sub_nodtyp(iend(iseg)) == 1)   {
                        if (nodpress(iend(iseg)) > nodpress(ista(iseg)))    {
                            disc_perf(x,y,z) += qq(iseg);
                            cntr += 1;
                        }
                    }
                }
                
                for (int inod = 0; inod < nnod; inod++) {
                    if (sub_nodtyp(inod) != nodtyp(inod) && sub_nodtyp(inod) != 0)   {
                        for (int iseg = 0; iseg < nseg; iseg++) {
                            if (ista(iseg) == inod && flagSeg(iseg) == 0)   {
                                if (nodpress(ista(iseg)) < nodpress(iend(iseg)))    {
                                    disc_perf(x,y,z) += qq(iseg);
                                    cntr += 1;
                                }
                            }
                            else if (iend(iseg) == inod && flagSeg(iseg) == 0)   {
                                if (nodpress(iend(iseg)) < nodpress(ista(iseg)))    {
                                    disc_perf(x,y,z) += qq(iseg);
                                    cntr += 1;
                                }
                            }

                        }
                    }
                }
                //disc_perf(x,y,z) /= cntr;
                
                // Boundary condition density
                for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                    if (flagNod(bcnod(inodbc)) == 1)    {
                        BC_dens(x,y,z) += 1;
                    }
                }


                double blockVol = sum(lseg(find(flagSeg == 1)) % (rseg(find(flagSeg == 1)) % rseg(find(flagSeg == 1))))*M_PI;
                // Vascular density
                block_dens(x,y,z) = blockVol/((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
                
                m = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*rho;
                disc_perf(x,y,z) = disc_perf(x,y,z)*(100/m)*1e-6;
            }
        }
    }
    
    FILE *ofp;
    
    string rootname = "Discrete_Flow/Sub-Block Density.txt";
    rootname = root + rootname;
    ofp = fopen(rootname.c_str(),"w");
    for (int x = 0; x < disc_x; x++) {
        for (int y = 0; y < disc_y; y++) {
            for (int z = 0; z < disc_z; z++) {
                fprintf(ofp,"%.5lf\t",block_dens(x,y,z));
            }
            fprintf(ofp,"\n");
        }
        fprintf(ofp,"\n\n\n");
    }
    
    fclose(ofp);
    
    
    FILE *ofp2;
    
    BC_dens /= nnodbc;
    rootname = "Discrete_Flow/BC Density.txt";
    rootname = root + rootname;
    ofp2 = fopen(rootname.c_str(),"w");
    for (int x = 0; x < disc_x; x++) {
        for (int y = 0; y < disc_y; y++) {
            for (int z = 0; z < disc_z; z++) {
                fprintf(ofp2,"%.5lf\t",BC_dens(x,y,z));
            }
            fprintf(ofp2,"\n");
        }
        fprintf(ofp2,"\n\n\n");
    }
    
    fclose(ofp2);

    

    FILE *ofperf;
    
    rootname = "Discrete_Flow/Sub-Block Perfusion.txt";
    rootname = root + rootname;
    ofperf = fopen(rootname.c_str(),"w");
    for (int x = 0; x < disc_x; x++) {
        for (int y = 0; y < disc_y; y++) {
            for (int z = 0; z < disc_z; z++) {
                fprintf(ofperf,"%lf\t",disc_perf(x,y,z));
            }
            fprintf(ofperf,"\n");
        }
        fprintf(ofperf,"\n\n\n");
    }
    
    fclose(ofperf);
    
    
    return disc_perf;
}
