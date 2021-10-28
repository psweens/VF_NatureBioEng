//
//  OutputNet.cpp
//  Vascular Flows 2.0
//
//  Outputs network file
//
//  Created by Paul Sweeney on 02/07/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

extern string netName,root;
extern int mxx,myy,mzz,nodsegm;
extern float alx,aly,alz,lb,maxl;

extern vec nodpress;

void outputNet(const string &filename, const int &nseg, const int &nnod, const int &nnodbc, const uvec &segname, const uvec &vesstyp, const umat &segnodname, const vec &diam, const vec &q, const vec &hd, const uvec &nodname, const mat &cnode, const uvec &bcnodname, const uvec &bctyp, const vec &bcprfl, const vec &bchd)    {

    FILE *ofp;
    
    string rootname = root + filename;
    ofp = fopen(rootname.c_str(),"w");
    fprintf(ofp,"%s",netName.c_str());
    fprintf(ofp,"%f %f %f  Box dimensions in microns \n",alx,aly,alz);
    fprintf(ofp,"%i %i %i  No. of tissue points in x,y,z directions \n",mxx,myy,mzz);
    fprintf(ofp,"%f    Outer bound distance \n",lb);
    fprintf(ofp,"%f    Max. segment length \n",maxl);
    fprintf(ofp,"%i    Max. segments per node \n",nodsegm);
    fprintf(ofp,"%i    Total number of segments\n",nseg);
    fprintf(ofp,"Segname Type Start End Diam Flow[qL/min] Hd\n");
    for(int iseg = 0; iseg < nseg; iseg++)     {
        fprintf(ofp,"%lli %lli %lli %lli %f %f %f\n",segname(iseg),vesstyp(iseg),segnodname(0,iseg),segnodname(1,iseg),diam(iseg),q(iseg),hd(iseg));
    }
    fprintf(ofp,"%i Total Number of nodes\n", nnod);
    fprintf(ofp,"Nodname x y z\n");
    for(int inod = 0; inod < nnod; inod++)  {
        fprintf(ofp, "%lli %f %f %f %f\n",nodname(inod),cnode(0,inod),cnode(1,inod),cnode(2,inod),nodpress(inod));
    }
    fprintf(ofp,"%i Total Number of boundary nodes\n", nnodbc);
    fprintf(ofp,"Bcnodname Bctyp Bcprfl BcHd\n");
    for (int inodbc = 0; inodbc < nnodbc; inodbc++)    {
        fprintf(ofp,"%lli %lli %f %f\n",bcnodname(inodbc),bctyp(inodbc),bcprfl(inodbc),bchd(inodbc));
    }
    
    fclose(ofp);

}
