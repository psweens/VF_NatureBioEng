//
//  flow_menu.cpp
//  Vascular Flows 2.0
//
//  Flow menu interface
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <chrono>
#include <random>
#include <ctime>

#include "global_variables.h"
#include "global_params.h"
#include "discrete_flow.h"
#include "output_data.h"
#include "initial_setup.h"

using namespace std;

void flow_menu()    {
    
rewind1:;
    cout<<"\t\tRun discrete flow solver? (y/n) " ;
    string input1;
    cin>>input1;
    if (input1 == "y")  {
        
    rewind2:;
        cout<<"\t\tRandomly assign flow boundary conditions? (y/n) ";
        string input2;
        cin>>input2;
        if (input2 == "y")  {
            cout<<"\t\tNo. of boundary conditions to apply? ";
            int input3;
            cin>>input3;
            assignBCs(input3);
        }
        else {
            cout<<"\t\t... No flow conditions applied."<<endl;
        }
        
        clock_t run_start;
        run_start = clock();
        
        
        // Setting up flow arrays
        flowArrays();
        
        
        // Allows for multiple runs of flow solver if BCs are assigned randomly
        int num = 12;
        string save_root = root;
        for (int i = 0; i < num; i++)    {
            
            q = zeros<vec>(nseg);
            flow_script(num);
            
            root = save_root;
        }
        
        
        
    }
    else if (input1 == "n") {
    rewind3:;
        cout<<"\t\tUse flow data from network file? (y/n) " ;
        string input4;
        cin>>input4;
    }
    else if (input1 != "n") {
        cout<<"\t\t*** Error: Input not within predefined parameters ***"<<endl;
        goto rewind1;
    }
    
}


void flow_script(const int &num)  {
    
    
    // Pressure BCs
    if (num > 1)    {
        randomBC(num);
    }
    
    
    // Outputs pre-solved network file
    outputNet("Discrete_Flow/Discrete_Output.txt", nseg, nnod, nnodbc, segname, vesstyp, segnodname, diam, q, hd, nodname, cnode, bcnodname, bctyp, bcprfl, bchd);
    
    storeBC = bcprfl;
    storeBCtyp = bctyp;
    storeBChd  = bchd;
    vec storeHD = hd;
    
    if (any(bctyp == 3))    {
        cout<<"\t\t\t*** Unknown boundary condition(s) detected ***"<<endl;
    }
    
    
    if (phaseseparation == 1)   {
        cout<<"\t\t\t...Incorporating variable hematocrit"<<endl;
    }
    
    
    flowSetup();
    cout<<"\t\t\tCalculating shear stress..."<<endl;
    int check = 0;
    int flowest = 0;
    // Flow estimation loop
    vec storeZero = zeros<vec>(nseg);
    while (check == 0)    {
        
        if (any(bctyp == 3))    {
            flowest = 1;
        }
        else    {
            check = 1;
        }
        
        oldFlowSign = flowSign;
        oldTau = tau;
        
        bcprfl = storeBC;
        bctyp = storeBCtyp;
        bchd = storeBChd;
        hd = storeHD;
        
        RBF = 0.;
        inflow = 0.;
        
        
        flow(nseg, nnod, nnodbc, segname, nodname, bcnodname, cnode, ista, iend, diam, lseg, q, qold, hd, hdold, nodtyp, bctyp, bcnod, bchd, bcprfl, c, conduc, segpress, nodpress, qq, tau, BCpress, BCflow);
        
        
        flowSign = sign(q);
        if (mean(abs(tau)) > 0.0)   {
            tau0 = mean(abs(tau))*beta*flowSign;
        }
        else    {
            tau0 %= flowSign;
        }
        // Compensating for zero flows - Fry is unlikely to converge otherwise
        tau0(find(q == 0.0)).fill(0.0);
        oldFlowSign(find(q == 0.0)).fill(0.0);
        
        
        int counter = 0;
        storeZero(find(zeroFlowVec == 1)).fill(1);
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (flowSign(iseg) != oldFlowSign(iseg) && storeZero(iseg) == 0)    {
                counter += 1;
            }
        }
        if (counter > 0 && counter <= sum(zeroFlowVec)){
            counter = 0;
        }

        
        oldktau = ktau;
        if (counter == 0 && ktau > 1e3 && flowest == 1) {
            check = 1;
            RBF = 0.;
            inflow = 0.;
            ktau /= 2;
            targStress = mean(abs(oldTau));
            tau0 = targStress*beta*sign(oldFlowSign);
            tau0(find(q == 0.0)).fill(0.0);
            flow(nseg, nnod, nnodbc, segname, nodname, bcnodname, cnode, ista, iend, diam, lseg, q, qold, hd, hdold, nodtyp, bctyp, bcnod, bchd, bcprfl, c, conduc, segpress, nodpress, qq, tau, BCpress, BCflow);
            cout<<"\n\t\t\t\tFinal ktau = "<<ktau<<"\tTarget wall shear stress = "<<targStress<<" dyn/cm2"<<endl;
        }
        
        
        ktau *= 2;
        
    }
    

    // Calculate perfusion
    vasc_perf(2);
    
    
    // Calculate segment pressures
    for (int iseg = 0; iseg < nseg; iseg++) {
        segpress(iseg) = (nodpress(ista(iseg)) + nodpress(iend(iseg)))/2.;
    }
    
    
    // Store pressures at boundary nodes
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        BCpress(inodbc) = nodpress(bcnod(inodbc));
    }
    
    
    run_time = (clock()-run_start)/(double) CLOCKS_PER_SEC;
    
    if (run_time > 60)    {
        cout<<"\t\t\t\tDiscrete flow solver run time: "<<run_time/60<<" minutes"<<endl;
    }
    else    {
        cout<<"\t\t\t\tDiscrete flow solver run time: "<<run_time<<" seconds"<<endl;
        sec_min = 1;
    }
    
    
    // Output data
    flowOutput();

    
    
}
