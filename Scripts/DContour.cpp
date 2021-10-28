//
//  DContour.cpp
//  Vascular Flows 2.0
//
//  Taken from Greens O2 by Tim Secomb
//  and update for the armadillo library
//
//  Created by Paul Sweeney on 16/07/2015.
//

#include <stdio.h>
#include <stdlib.h>

#include "global_variables.h"
#include "global_params.h"
#include "inter_func.h"
#include "inter_var.h"
#include "misc_func.h"

void DContour(const string &fname, const int &flag)   {
    
    string filename = root + fname;

    int isl1,isl2;
    double xs,ys,xmax,ymax,xmin,ymin;
    double red,green,blue,pressureplot = 0.0;
    
    FILE *ofp;
    
    xmin = min(cnode.row(0));
    xmax = max(cnode.row(0));
    ymin = min(cnode.row(1));
    ymax = max(cnode.row(1));
    
    picfac = 500./fmax(xmax,ymax);


    
    xsl0(0) = xmin;
    xsl0(1) = ymin;
    xsl1(0) = xmax;
    xsl1(1) = ymin;
    xsl2(0) = xmin;
    xsl2(1) = ymax;
    
    
    int nsl1 = (alx/100);
    int nsl2 = (aly/100);
    
    int nlmax = 1;
    int NL = 60;
    double pint = 0.0;
    if (flag == 0)  {
        if (iP < segpress.min())    {
            pint = (segpress.max()-iP)/NL;
        }
        else {
            pint = (segpress.max()-segpress.min())/NL;
        }
    }

    
    if(NL > nlmax) nlmax = NL;
    
    vec cl = zeros<vec>(nlmax+1);
    mat zv = zeros<mat>(nsl1+1,nsl2+1);
    
    vec x = zeros<vec>(3);
    
    printf("\t\t...Generating contour plots\n");
    //Calculate P on a planar slice through the region, specified by three corners and number of points along each edge
    int cntr = 0;
    for(isl1 = 0; isl1 <= nsl1; isl1++) {
        for(isl2 = 0; isl2 <= nsl2; isl2++){
            for(int i = 0; i < 3; i++)  {
                x(i) = xsl0(i) + isl1*(xsl1(i)-xsl0(i))/nsl1 + isl2*(xsl2(i)-xsl0(i))/nsl2;
            }
            double p = 0.0;
            if (flag == 0)  {
                p = DCeval(x,flag) + iP;
            }
            else if (flag == 1) {
                p = log(abs(DCeval(x,flag))*1e3); // log(um / s)
            }
            cntr += 1;
            //printf("p=%f\n",p);
            zv(isl1,isl2) = p;
            //printf("zv(%i,%i)=%f\n",isl1,isl2,zv(isl1,isl2));
        }
    }
    double vmin = 0.0;
    double vmax = 0.0;
    if (flag == 1)  {
        vmin = zv.min();
        vmax = zv.max();
        pint = (vmax - vmin) / NL;
        
        outputf(ift,"Min. Interstitial Velocity", exp(vmin),"um/s");
        outputf(ift,"Max. Interstitial Velocity", exp(vmax),"um/s");
    }
    else {
        if (iP < segpress.min())    {
            vmin = iP;
        }
        else {
            vmin = segpress.min();
        }
        vmax = segpress.max();
        vmin = zv.min();
        vmax = zv.max();
        pint = (vmax - vmin) / NL;
        outputf(ift,"Min. Interstitial Pressure", zv.min(),"mmHg");
        outputf(ift,"Max. Interstitial Pressure", zv.max(),"mmHg");
    }

    
    ofp = fopen(filename.c_str(), "w");
    fprintf(ofp, "%%!PS-Adobe-2.0\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "12 scalefont\n");
    fprintf(ofp, "setfont\n");
    for(int i = 0; i <= NL; i++) cl(i) = vmin + (i)*pint;
    

    DContour_shade(ofp,nsl1,nsl2,picfac,NL,xmin,xmax,ymin,ymax,cl,zv);

    
    // Plot bounding box
    
    fprintf(ofp, "/mx {%g mul 50 add} def\n",picfac);
    fprintf(ofp, "/my {%g mul 100 add} def\n",picfac);
    
    fprintf(ofp, "newpath\n");
    fprintf(ofp, "%g mx %g my m\n",xmin-xmin,ymin-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmax-xmin,ymin-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmax-xmin,ymax-ymin);
    fprintf(ofp, "%g mx %g my l\n",xmin-xmin,ymax-ymin);
    fprintf(ofp, "closepath\n");
    fprintf(ofp, "stroke\n");
    
    
    if (flag == 2)  {
        // Plots projection of the network colour coded according to average segment pressure
        //printf("Pa=%f\n",Pa);
        
        fprintf(ofp,"0 0 0 setrgbcolor\n");//black
        for(int iseg = 0; iseg < nseg; iseg++){
            fprintf(ofp,"%g setlinewidth\n",picfac*(diam(iseg)+20.)); // extra 20 on diameter bounds the vessels in black
            xs = cnode(0,ista(iseg)) -xmin;
            ys = cnode(1,ista(iseg)) -ymin ;
            fprintf(ofp, "%g mx %g my m ", xs,ys);
            xs = cnode(0,iend(iseg)) -xmin ;
            ys = cnode(1,iend(iseg)) -ymin ;
            fprintf(ofp, "%g mx %g my l ", xs,ys);
            fprintf(ofp, "stroke\n");
        }
        
        
        for(int iseg = 0; iseg < nseg; iseg++){
            if (flag == 0)     {
                pressureplot = (segpress(iseg) - segpress.min())/(segpress.max()-segpress.min());
                if (iP < segpress.min())    {
                    pressureplot = (segpress(iseg) - iP)/(segpress.max()-iP);
                }
                else {
                    pressureplot = (segpress(iseg) - segpress.min())/(segpress.max()-segpress.min());
                }
            }
            else if (flag == 1) {
                pressureplot = (vel(iseg) - vel.min())/(vel.max()-vel.min());
            }
            
            blue = fmin(fmax(1.5-4*fabs(pressureplot-0.25),0.),1.);
            green = fmin(fmax(1.5-4*fabs(pressureplot-0.5),0.),1.);
            red = fmin(fmax(1.5-4*fabs(pressureplot-0.75),0.),1.);
            if(pressureplot < 0. || pressureplot > 1.){
                red = 1.0;
                green = 1.0;
                blue = 1.0;
            }
            
            
            fprintf(ofp,"%6.3f %6.3f %6.3f setrgbcolor\n",red,green,blue);
            //		fprintf(ofp,"1 setlinewidth\n"); // vessels have uniform diameter
            fprintf(ofp,"%g setlinewidth\n",picfac*diam(iseg)); // vessels are drawn as proportional to their actual diameters
            xs = cnode(0,ista(iseg)) - xmin;
            ys = cnode(1,ista(iseg)) - ymin;
            fprintf(ofp, "%g mx %g my m ", xs,ys);
            xs = cnode(0,iend(iseg)) - xmin;
            ys = cnode(1,iend(iseg)) - ymin;
            fprintf(ofp, "%g mx %g my l ", xs,ys);
            fprintf(ofp, "stroke\n");
        }
    }
    
    
    // Create a colour bar
    double cbbox = 21.;    //size of boxes
    double cbx = 550;  //origin of color bar
    double cby = 50;   //origin of color bar
    vec bluevect = zeros<vec>(NL/3+1);
    vec greenvect = zeros<vec>(NL/3+1);
    vec redvect = zeros<vec>(NL/3+1);
    vec pvals = zeros<vec>(NL/3+1);
    for(int k = 0; k <=NL/3; k++){
        float xz = 3*float(k)/float(NL);
        bluevect(k) = min(max(1.5-4*abs(xz-0.25), 0.), 1.);
        greenvect(k)= min(max(1.5-4*abs(xz-0.5), 0.), 1.);
        redvect(k)  = min(max(1.5-4*abs(xz-0.75), 0.), 1.);
        pvals(k) = vmin + 3*k*(vmax - vmin) / NL;
    }
    fprintf(ofp,"%d setlinewidth\n",1);
    for(int k = 0; k <= NL/3; k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",redvect(k),greenvect(k),bluevect(k));
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
        fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%.3g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),pvals(k));
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(NL/3+1),cbx,cby+cbbox*(NL/3+1));

    
    
    fprintf(ofp, "showpage\n");
    fclose(ofp);
    printf("\t\t\t ...Done\n");
    
}
