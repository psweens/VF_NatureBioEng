//
//  DContour_shade.cpp
//  Vascular Flows 2.0
//
//  Taken from Greens O2 by Tim Secomb
//  and update for the armadillo library
//
//  Created by Paul Sweeney on 16/07/2015.
//

#include <stdio.h>
#include <armadillo>
#include <cmath>

#include "global_variables.h"
#include "inter_var.h"

using namespace arma;


void DContour_shade(FILE *ofp, const int &m, const int &n, double &scalefac, int &nl, const double &xmin, const double &xmax, const double &ymin, const double &ymax, const vec &cl, const mat &zv)  {
    
    fprintf(ofp, "%%!PS\n");
    int iwsp,in1,in2,in3,in4,inh;
    const int ncmax=1e6;
    double xz = 0.0,xv12 = 0.0,yv12 = 0.0,xv23 = 0.0,yv23 = 0.0,xv34 = 0.0,yv34 = 0.0,xv41 = 0.0,yv41 = 0.0;
    double cx,cy,cbx,cby,cbbox,eps;
    
    cx = 50;    //25; //origin of contour plot
    cy = 100;   //50;//origin of contour plot
    cbbox = 21.;    //size of boxes
    cbx = 550;  //origin of color bar
    cby = 50;   //origin of color bar
    
    // vectors were indexed from 0 to nl --> nl+1 entries so we now index from 0 to nl+1
    vec xv = zeros<vec>(m+1);
    vec yv = zeros<vec>(n+1);
    vec xvv = zeros<vec>(5);
    vec yvv = zeros<vec>(5);
    vec red = zeros<vec>(nl+1);
    vec green = zeros<vec>(nl+1);
    vec blue = zeros<vec>(nl+1);
    mat wsp = zeros<mat>(ncmax+1,5);
    imat corners = zeros<imat>(5,2);
    vec dd = zeros<vec>(5);
    ivec ii = zeros<ivec>(5);
    ivec jj = zeros<ivec>(5);
    corners(1,0) = 0;
    corners(1,1) = 0;
    corners(2,0) = 1;
    corners(2,1) = 0;
    corners(3,0) = 1;
    corners(3,1) = 1;
    corners(4,0) = 0;
    corners(4,1) = 1;
    
    
    for(int i=1; i<=m; i++) xv(i) = xmin + (i - 1)*(xmax - xmin)/(m -1) - xmin;
    for(int j=1; j<=n; j++) yv(j) = ymin + (j - 1)*(ymax - ymin)/(n -1) - ymin;
    
    fprintf(ofp, "/mx {%g sub %g mul %g add} def\n",xmin-xmin,scalefac,cx);
    fprintf(ofp, "/my {%g sub %g mul %g add} def\n",ymin-ymin,scalefac,cy);
    fprintf(ofp, "/m {moveto} def\n");
    fprintf(ofp, "/l {lineto} def\n");
    fprintf(ofp, "/n {newpath} def\n");
    fprintf(ofp, "/s {stroke} def\n");
    fprintf(ofp, "/cf {closepath fill} def\n");
    fprintf(ofp, "/cs {closepath stroke} def\n");
    fprintf(ofp, "0.5 setlinewidth\n");
    fprintf(ofp, "/Times-Roman findfont\n");
    fprintf(ofp, "12 scalefont\n");
    fprintf(ofp, "setfont\n");
    
    //eps defines extra width added to colored rectangles.
    //otherwise, small lines of different color appear between rectangles
    eps = fmax(xmax-xmin,ymax-ymin)*1.e-3;
    
    //Set up colors using Matlab 'jet' scheme
    for(int k=0; k<=nl; k++){
        xz = float(k)/float(nl);
        blue(k) = fmin(fmax(1.5-4*fabs(xz-0.25), 0.), 1.);
        green(k)= fmin(fmax(1.5-4*fabs(xz-0.5), 0.), 1.);
        red(k)  = fmin(fmax(1.5-4*fabs(xz-0.75), 0.), 1.);
    }
    //Color whole region with lowest color
    fprintf(ofp, "%f %f %f setrgbcolor\n",red(0),green(0),blue(0));
    fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cf\n",
            xmin-xmin,ymin-ymin,xmax-xmin,ymin-ymin,xmax-xmin,ymax-ymin,xmin-xmin,ymax-ymin);
    //Analyze each rectangle separately. Overwrite lower colors
    iwsp = 0;
    for(int k=0; k<=nl; k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",red(k),green(k),blue(k));
        for(int i=0; i<m; i++) for(int j=1; j<n; j++){
            in1 = 0;
            for(int in=0; in<=4; in++){
                ii(in) = i + corners(in,0);
                jj(in) = j + corners(in,1);
                dd(in) = zv(ii(in),jj(in)) - cl(k);
                if(dd(in) >= 0.) in1 = in;	//find a corner this color or higher
            }
            inh = 0;
            if(k < nl) for(int in=0; in<=4; in++)
                if(zv(ii(in),jj(in)) - cl(k+1) < 0.) inh = 1;//check that not all corners are above the next contour
            if(in1 > 0 && inh > 0){
                in2 = in1%4 + 1;
                in3 = in2%4 + 1;
                in4 = in3%4 + 1;
                for(int in=1; in<=4; in++){
                    xvv(in) = xv(ii(in));
                    yvv(in) = yv(jj(in));
                    if(ii(in) == i+1) xvv(in) += eps; else xvv(in) -= eps;
                    if(jj(in) == j+1) yvv(in) += eps; else yvv(in) -= eps;
                }
                if(dd(in1) != dd(in2)){
                    xv12 = (dd(in1)*xv(ii(in2))-dd(in2)*xv(ii(in1)))/(dd(in1)-dd(in2));
                    yv12 = (dd(in1)*yv(jj(in2))-dd(in2)*yv(jj(in1)))/(dd(in1)-dd(in2));
                }
                if(dd(in2) != dd(in3)){
                    xv23 = (dd(in2)*xv(ii(in3))-dd(in3)*xv(ii(in2)))/(dd(in2)-dd(in3));
                    yv23 = (dd(in2)*yv(jj(in3))-dd(in3)*yv(jj(in2)))/(dd(in2)-dd(in3));
                }
                if(dd(in3) != dd(in4)){
                    xv34 = (dd(in3)*xv(ii(in4))-dd(in4)*xv(ii(in3)))/(dd(in3)-dd(in4));
                    yv34 = (dd(in3)*yv(jj(in4))-dd(in4)*yv(jj(in3)))/(dd(in3)-dd(in4));
                }
                if(dd(in4) != dd(in1)){
                    xv41 = (dd(in4)*xv(ii(in1))-dd(in1)*xv(ii(in4)))/(dd(in4)-dd(in1));
                    yv41 = (dd(in4)*yv(jj(in1))-dd(in1)*yv(jj(in4)))/(dd(in4)-dd(in1));
                }
                fprintf(ofp, "n %g mx %g my m ",xvv(in1),yvv(in1));
                if(dd(in2) > 0){													//corners 1,2 are this color
                    fprintf(ofp, "%g mx %g my l ",xvv(in2),yvv(in2));
                    if(dd(in3) > 0){												//corners 1,2,3 are this color
                        fprintf(ofp, "%g mx %g my l ",xvv(in3),yvv(in3));
                        if(dd(in4) > 0)
                            fprintf(ofp, "%g mx %g my l ",xvv(in4),yvv(in4));		//corners 1,2,3,4 are this color
                        else{														//corners 1,2,3,not 4 are this color
                            fprintf(ofp, "%g mx %g my l ",xv34,yv34);
                            fprintf(ofp, "%g mx %g my l ",xv41,yv41);
                            iwsp += 1;
                            wsp(iwsp,1) = xv34;
                            wsp(iwsp,2) = yv34;
                            wsp(iwsp,3) = xv41;
                            wsp(iwsp,4) = yv41;
                        }
                    }
                    else{															//corners 1,2,not 3 are this color
                        fprintf(ofp, "%g mx %g my l ",xv23,yv23);
                        iwsp += 1;
                        wsp(iwsp,1) = xv23;
                        wsp(iwsp,2) = yv23;
                        if(dd(in4) > 0){											//corners 1,2,not 3,4 are this color
                            fprintf(ofp, "%g mx %g my l ",xv34,yv34);
                            wsp(iwsp,3) = xv34;
                            wsp(iwsp,4) = yv34;
                            fprintf(ofp, "%g mx %g my l ",xvv(in4),yvv(in4));
                        }
                        else{
                            fprintf(ofp, "%g mx %g my l ",xv41,yv41);				//corners 1,2,not 3,not 4 are this color
                            wsp(iwsp,3) = xv41;
                            wsp(iwsp,4) = yv41;
                        }
                    }
                }
                else{																//corners 1,not 2 are this color
                    fprintf(ofp, "%g mx %g my l ",xv12,yv12);
                    iwsp += 1;
                    wsp(iwsp,1) = xv12;
                    wsp(iwsp,2) = yv12;
                    if(dd(in3) > 0){												//corners 1,not 2,3 are this color
                        fprintf(ofp, "%g mx %g my l ",xv23,yv23);
                        wsp(iwsp,3) = xv23;
                        wsp(iwsp,4) = yv23;
                        fprintf(ofp, "%g mx %g my l ",xvv(in3),yvv(in3));
                        if(dd(in4) > 0)
                            fprintf(ofp, "%g mx %g my l ",xvv(in4),yvv(in4));		//corners 1,not 2,3,4 are this color
                        else{														//corners 1,not 2,3,not 4 are this color
                            fprintf(ofp, "%g mx %g my l ",xv34,yv34);
                            fprintf(ofp, "%g mx %g my l ",xv41,yv41);
                            iwsp += 1;
                            wsp(iwsp,1) = xv34;
                            wsp(iwsp,2) = yv34;
                            wsp(iwsp,3) = xv41;
                            wsp(iwsp,4) = yv41;
                        }
                    }
                    else{															//corners 1,not 2,not 3 are this color
                        if(dd(in4) > 0){											//corners 1,not 2,not 3,4 are this color
                            fprintf(ofp, "%g mx %g my l ",xv34,yv34);
                            wsp(iwsp,3) = xv34;
                            wsp(iwsp,4) = yv34;
                            fprintf(ofp, "%g mx %g my l ",xvv(in4),yvv(in4));
                        }
                        else{
                            fprintf(ofp, "%g mx %g my l ",xv41,yv41);				//corners 1,not 2,not 3,not 4 are this color
                            wsp(iwsp,3) = xv41;
                            wsp(iwsp,4) = yv41;
                        }
                    }
                }
                if(iwsp > ncmax-4) printf("*** Error: ncmax too small in contr\n");
                fprintf(ofp, "cf\n");
            }
        }
    }
    //Now outline contours
    /*fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    for(int in=1; in<=iwsp; in++) fprintf(ofp, "n %g mx %g my m %g mx %g my l s\n",
                                      wsp(in,1),wsp(in,2),wsp(in,3),wsp(in,4));*/
    //Draw a box
    fprintf(ofp, "0 0 0 setrgbcolor\n");//black
    fprintf(ofp, "n %g mx %g my m %g mx %g my l %g mx %g my l %g mx %g my l cs\n",
            xmin-xmin,ymin-ymin,xmax-xmin,ymin-ymin,xmax-xmin,ymax-ymin,xmin-xmin,ymax-ymin);
    

    
    /*fprintf(ofp,"%d setlinewidth\n",1);
    for(int k=0; k<=nl; k++){
        fprintf(ofp, "%f %f %f setrgbcolor\n",red(k),green(k),blue(k));
        fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cf\n",
                cbx,cby+k*cbbox,cbx+cbbox,cby+k*cbbox,cbx+cbbox,cby+(k+1)*cbbox,cbx,cby+(k+1)*cbbox);
        //if(k>0) {
            fprintf(ofp, "%g %f m 0 0 0 setrgbcolor (%.3g) show\n",cbx+cbbox*1.1,cby+cbbox*(k-0.1),cl(k));
            cout<<cl(k)<<endl;
        //}
    
    }
    fprintf(ofp, "n %g %g m %g %g l %g %g l %g %g l cs\n",
            cbx,cby,cbx+cbbox,cby,cbx+cbbox,cby+cbbox*(nl+1),cbx,cby+cbbox*(nl+1));*/
}
