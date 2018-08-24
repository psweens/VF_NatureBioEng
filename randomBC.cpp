//
//  randomBC.cpp
//  Vascular-Flow
//
//  Randomly assigns pressure boundary conditions
//  Reads in convex hull of network using 'hull_read' function
//
//  Created by Paul Sweeney on 19/04/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>

#include "global_params.h"
#include "global_variables.h"
#include "output_data.h"

void randomBC(const int &check) {
    
    
    // Assigning pressure boundary conditions
    string bcInput;
    if (check > 1)  {
        bcInput = "y";
        goto here;
    }
    cout<<"\n\t\t\tRandomly assign pressure boundary conditions? ";
    cin>>bcInput;
here:;
    if (bcInput == "y") {

        unsigned seed = (double) chrono::system_clock::now().time_since_epoch().count();
        ivec randBC = zeros<ivec>(nnodbc);
        for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
            randBC(inodbc) = inodbc;
        }
        
        shuffle(&randBC(0), &randBC(randBC.n_elem-1), default_random_engine(seed)); // Randomly shuffle boundary node indexes
        bcprfl.fill(0.);
        bctyp.fill(3);
        
        int BCran = 0;
        srand((double) time(0));
        float Acntr = 0.;
        float Vcntr = 0.;
        
        // Locates BC nodes based on convex hull - assigns values based on exclusion region around each node
        ivec flag = zeros<ivec>(nnodbc);
        mat temp;
        hull_read("Read_Hull",flag,temp,0);
        int i = 0;
        int j = 0;
        high = 30;
        low = 20;
        double exclude = 0.05;
        while (i < round(0.05*nnodbc))  {
            if (flag(randBC(j)) == 1)   {
                bctyp(randBC(j)) = 0;
                BCran = rand()%2;
                if (BCran == 1)   {
                    int run_away = 0;
                    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                        if (inodbc != randBC(j) && abs(cnode(0,bcnod(inodbc)) - cnode(0,bcnod(randBC(j)))) < exclude*alx && abs(cnode(1,bcnod(inodbc)) - cnode(1,bcnod(randBC(j)))) < exclude*aly && abs(cnode(2,bcnod(inodbc)) - cnode(2,bcnod(randBC(j)))) < exclude*alz && bcprfl(inodbc) == low) {
                            run_away = 1;
                        }
                    }
                    if (run_away == 1)  {
                        bcprfl(randBC(j)) = low;
                        Vcntr += 1;
                    }
                    else {
                        bcprfl(randBC(j)) = high;
                        Acntr += 1;
                    }
                }
                else {
                    int run_away = 0;
                    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
                        if (inodbc != randBC(j) && abs(cnode(0,bcnod(inodbc)) - cnode(0,bcnod(randBC(j)))) < exclude*alx && abs(cnode(1,bcnod(inodbc)) - cnode(1,bcnod(randBC(j)))) < exclude*aly && abs(cnode(2,bcnod(inodbc)) - cnode(2,bcnod(randBC(j)))) < exclude*alz && bcprfl(inodbc) == high) {
                            run_away = 1;
                        }
                    }
                    if (run_away == 1)  {
                        bcprfl(randBC(j)) = high;
                        Acntr += 1;
                    }
                    else {
                        bcprfl(randBC(j)) = low;
                        Vcntr += 1;
                    }
                }
                i += 1;
            }
            j += 1;
        }
        i = 0;
        j = 0;
        double no_flow = 0.33;
        while (i < round(no_flow*nnodbc))  {
            if (flag(randBC(j)) == 0)   {
                bctyp(randBC(j)) = 1;
                bcprfl(randBC(j)) = 0.0;
                i += 1;
            }
            j += 1;
        }
        
        perc1 = Acntr/nnodbc*100;
        perc2 = Vcntr/nnodbc*100;
        printf("\t\t\t%.0f mmHg assigned to %.1f%% of BCs\n",high,perc1);
        printf("\t\t\t%.0f mmHG assigned to %.1f%% of BCs\n",low,perc2);
        
    }
    
}
