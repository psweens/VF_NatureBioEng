//
//  assignBCs.cpp
//  Vascular Flows 2.0
//
//  Randomly assigning boundary conditions based on data from Simon, Tim Secomb & David Boas
//  flow solutions
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <chrono>
#include <random>
#include <ctime>
#include <cmath>

#include "discrete_flow.h"
#include "global_variables.h"

using namespace std;

void assignBCs(const int &numBC)    {
    
    dataFitFlag = zeros<ivec>(nnodbc);
    
    unsigned seed = (double) chrono::system_clock::now().time_since_epoch().count();
    
    int nFlowPress = 0;
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 0 || bctyp(inodbc) == 1)   {
            nFlowPress += 1;
        }
    }

    ivec randBC = zeros<ivec>(nnodbc-nFlowPress);
    
    int i = 0;
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 3) {
            randBC(i) = inodbc;
            i += 1;
        }
    }
    
    shuffle(&randBC(0), &randBC(randBC.n_elem-1), default_random_engine(seed)); // Randomly shuffle boundary node indexes

    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (bcnod(inodbc) == ista(iseg))    {
                for (int jseg = 0; jseg < nseg; jseg++ )    {
                    if (iend(iseg) == ista(jseg) || iend(iseg) == iend(jseg))   {
                        //cout<<diam(iseg)<<"\t"<<diam(jseg)<<endl;
                    }
                }
            }
            else if (bcnod(inodbc) == iend(iseg))    {
                for (int jseg = 0; jseg < nseg; jseg++ )    {
                    if (ista(iseg) == ista(jseg) || ista(iseg) == iend(jseg))   {
                        //cout<<diam(iseg)<<"\t"<<diam(jseg)<<endl;
                    }
                }
            }
        }
    }
    int cntr = 0;
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (diam(iseg) < 110. && (bcnod(inodbc) == ista(iseg) || bcnod(inodbc) == iend(iseg)) && bctyp(inodbc) == 3) {
                cntr += 1;
            }
        }
    }
    
    i = 0;
    cntr = 0;
    while (i < nnodbc-nFlowPress)   {
        for (int iseg = 0; iseg < nseg; iseg++)    {
            if (diam(iseg) < 110. && (bcnod(randBC(i)) == ista(iseg) || bcnod(randBC(i)) == iend(iseg)))  {
                int ran = rand()%2;
                if (ran == 0)   {
                    ran = -1;
                }
                // Data-fitted equation for Tim's data
                //bcprfl(inodbc) = 0.1351*exp(0.4325*diam(iseg))*ran;
                // Combined Secomb/Walker-Samuel/Boas datasets
                bcprfl(randBC(i)) = 0.0513*powf(diam(iseg),1.937)*ran;
                bctyp(randBC(i)) = 1;
                dataFitFlag(randBC(i)) = 1;
                cntr += 1;
            }
        }
        if (cntr == numBC)  {
            goto exit;
        }
        i += 1;
    }
exit:;

    int calc = cntr*100/nnodbc;
    if (cntr != numBC)  {
        cout<<"\t\t\t***Warning: Flow conditions could only be applied to "<<calc<<"% of unknown boundary conditions ***"<<endl;
    }
    
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 0) {
            dataFitFlag(inodbc) = 2;
        }
    }
    

    ivec geoBC = zeros<ivec>(nnodbc);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (bcnod(inodbc) == ista(iseg) || bcnod(inodbc) == iend(iseg)) {
                if (vesstyp(iseg) == 1)    {
                    geoBC(inodbc) = 1;
                }
                else if (vesstyp(iseg) == 2)   {
                    geoBC(inodbc) = 2;
                }
                else if (vesstyp(iseg) == 3)   {
                    geoBC(inodbc) = 3;
                }
            }
        }
    }
    double art = 0;
    double cap = 0;
    double ven = 0;
    int cnt = 0;
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        if (bctyp(inodbc) == 1) {
            cnt += 1;
            if (geoBC(inodbc) == 1)  {
                art += 1;
            }
            else if (geoBC(inodbc) == 2)    {
                cap += 1;
            }
            else if (geoBC(inodbc) == 3)    {
                ven += 1;
            }
        }
    }
    //cout<<cnt<<endl;
    art /= numBC;
    cap /= numBC;
    ven /= numBC;
    art *= 100;
    cap *= 100;
    ven *= 100;
    double perc = art + cap + ven;
    
    /*cout<<"\t\tPercentage of assigned BCs: "<<perc<<"%"<<endl;
    cout<<"\t\t\tPer Geometry: "<<endl;
    cout<<"\t\t\t\tArteriole: "<<art<<"%"<<endl;
    cout<<"\t\t\t\tCapillary: "<<cap<<"%"<<endl;
    cout<<"\t\t\t\tVenule: "<<ven<<"%"<<endl;*/
    
}
