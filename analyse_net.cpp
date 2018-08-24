//
//  analyseNet.cpp
//  Vascular Flows 2.0
//
//  Created by Paul Sweeney on 03/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <cmath>

#include "global_variables.h"
#include "global_params.h"

using namespace arma;
using namespace std;

void analyse_net()   {
    
    printf("\t No. of segments = %i\n",nseg);
    printf("\t No. of nodes = %i\n",nnod);
    printf("\t No. of boundary nodes = %i\n",nnodbc);
    
    for (int iseg = 0; iseg < nseg; iseg++)	{
        //Search for nodes corresponding to this segment
        for (int i = 0; i < 2; i++) {
            for (int inod = 0; inod < nnod; inod++) {
                if (nodname(inod) == segnodname(i,iseg))    {
                    if (i == 0) {
                        ista(iseg) = inod;
                        goto foundit;
                    }
                    else if (i == 1)    {
                        iend(iseg) = inod;
                        goto foundit;
                    }
                }
            }
            printf("*** Error: No matching node found for segname %lli\n", segname(iseg));
            foundit:;
        }
    }
    

    //Setup nodtyp, nodseg and nodnod
    for (int iseg = 0; iseg < nseg; iseg++) {
        int inod1 = (int) ista(iseg);
        int inod2 = (int) iend(iseg);
        nodtyp(inod1) += 1;
        nodtyp(inod2) += 1;
        if(nodtyp(inod1) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod1);
        }
        if(nodtyp(inod2) > nodsegm) {
            printf("*** Error: Too many segments connected to node %i\n", inod2);
        }
        /*nodseg(nodtyp(inod1),inod1) = iseg;
        nodseg(nodtyp(inod2),inod2) = iseg;
        nodnod(nodtyp(inod1),inod1) = inod2;
        nodnod(nodtyp(inod2),inod2) = inod1;*/
    }

    
    for (int inodbc = 0; inodbc < nnodbc; inodbc++){
        // Search for node corresponding to this node name
        for(int inod = 0; inod < nnod; inod++) {
            if(nodname(inod) == bcnodname(inodbc))  {
                bcnod(inodbc) = inod;
                if(nodtyp(inod) != 1)   {
                    printf("*** Error: Boundary node %lli is not a 1-segment node\n", nodname(inod));
                }
            goto foundit2;
            }
        }
        printf("*** Error: No matching node found for nodname %lli, inodbc %i\n", bcnodname(inodbc),inodbc);
        foundit2:;
    }
    
    
    // start(k,iseg) = coordinates of starting point of segment iseg
    // End(k,iseg) = coordinates of ending point of segment iseg
    tlength = 0.;
    compLseg = 1;
    for(int iseg = 0; iseg < nseg; iseg++) {
        rseg(iseg) = diam(iseg) / 2.0;
        qq(iseg) = abs(q(iseg));
        for(int k = 0; k < 3; k++){
            start(k,iseg) = cnode(k,ista(iseg));
            End(k,iseg) = cnode(k,iend(iseg));
            ss(k) = End(k,iseg) - start(k,iseg);
        }
        if (compLseg == 1)  {
            lseg(iseg) = sqrt(pow(ss(0),2) + pow(ss(1),2) + pow(ss(2),2));
        }
        for(int k = 0; k < 3; k++)   {
            scos(k,iseg) = ss(k) / lseg(iseg);
        }
        tlength += lseg(iseg);
    }
    
    BCgeo = zeros<uvec>(nnodbc);
    for (int inodbc = 0; inodbc < nnodbc; inodbc++) {
        for (int iseg = 0; iseg < nseg; iseg++) {
            if (bcnod(inodbc) == ista(iseg) || bcnod(inodbc) == iend(iseg)) {
                if (vesstyp(iseg) == 1)    {
                    BCgeo(inodbc) = 1;
                }
                else if (vesstyp(iseg) == 2)   {
                    BCgeo(inodbc) = 2;
                }
                else {
                    BCgeo(inodbc) = 3;
                }
            }
        }
    }
}
