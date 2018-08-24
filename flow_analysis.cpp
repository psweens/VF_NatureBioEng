//
//  flow_analysis.cpp
//  Vascular-Flow
//
//  Spits out network flow data
//
//  Created by Paul Sweeney on 01/07/2016.
//  Copyright © 2016 Paul Sweeney - University College London. All rights reserved.
//

#include <stdio.h>
#include "global_variables.h"
#include "global_params.h"
#include "output_data.h"
#include "discrete_flow.h"

void flow_analysis(const string &filename)    {
    
    cout<<"\t\t\t\tRegional blood flow: "<<RBF<<" ml/min/100g"<<endl;
    cout<<"\t\t\t\tNetwork Inflow: "<<inflow<<" nl/min"<<endl;
    
    vel = (qq*1e3/60)/(M_PI*pow(0.5*diam,2));
    
    vec Re = zeros<vec>(nseg);
    for (int iseg = 0; iseg < nseg; iseg++) {
        double visc = viscor((diam(iseg)*1e3),hd(iseg))*xi;
        Re(iseg) = rho*vel(iseg)*lseg(iseg)*1e-3 / visc;
    }
    
    cout<<"\t\t\t\tAll Vessels: "<<mean(vel)<<"\t"<<stddev(vel)<<" mm/s"<<endl;
    if (max(vesstyp) == 3)  {
        cout<<"\t\t\t\tArterioles: "<<mean(vel(find(vesstyp == 1)))<<"\t"<<stddev(vel(find(vesstyp == 1)))<<" mm/s"<<endl;
        cout<<"\t\t\t\tCapillaries: "<<mean(vel(find(vesstyp == 2)))<<"\t"<<stddev(vel(find(vesstyp == 2)))<<" mm/s"<<endl;
        cout<<"\t\t\t\tVenules: "<<mean(vel(find(vesstyp == 3)))<<"\t"<<stddev(vel(find(vesstyp == 3)))<<" mm/s"<<endl;
    }
    cout<<"\t\t\t\tktau/kp : "<<(ktau/2)/kp<<endl;
    
    vec vstt = zeros<vec>(nseg);
    vstt = 60*1e-6*(M_PI*pow(diam(find(qq > 0.0)),2) % lseg(find(qq > 0.0))) / (4*qq(find(qq > 0.0)));
    
    FILE *ofp;
    
    string rootname = root + filename;
    
    ofp = fopen(rootname.c_str(),"w");
    if (sec_min == 0)   {
        fprintf(ofp,"Algorithm Run-Time = %f sec\n\n",run_time);
    }
    else {
        fprintf(ofp,"Algorithm Run-Time = %f min\n\n",run_time);
    }
    fprintf(ofp,"%.1f mmHg : %.1f%%\n",high,perc1);
    fprintf(ofp,"%.1f mmHg : %.1f%%\n",low,perc2);
    fprintf(ofp,"Regional Blood Flow = %.3f ml/min/100g\n",RBF);
    fprintf(ofp,"Network Inflow = %.3f nl/min\n",inflow);
    fprintf(ofp,"Reynolds No.: Min. = %.3e, Max. =  %.3e, Mean = %.3e ± %.3e\n",min(Re),max(Re),mean(Re),stddev(Re));
    fprintf(ofp,"Mean Network Pressure = %.3f ± %.3f mmHg\n",mean(segpress),stddev(segpress));
    if (vesstyp.max() == 3) {
        fprintf(ofp,"\tMean Arteriolar Pressure = %.3f ± %.3f mmHg\n",mean(segpress(find(vesstyp == 1))),stddev(segpress(find(vesstyp == 1))));
        fprintf(ofp,"\tMean Capillary Pressure = %.3f ± %.3f mmHg\n",mean(segpress(find(vesstyp == 2))),stddev(segpress(find(vesstyp == 2))));
        fprintf(ofp,"\tMean Venular Pressure = %.3f ± %.3f mmHg\n",mean(segpress(find(vesstyp == 3))),stddev(segpress(find(vesstyp == 3))));
    }
    fprintf(ofp,"Mean Network Flow = %.3f ± %.3f nl/min\n",mean(qq),stddev(qq));
    if (vesstyp.max() == 3) {
        fprintf(ofp,"\tMean Arteriolar Flow = %.3f ± %.3f nl/min\n",mean(qq(find(vesstyp == 1))),stddev(qq(find(vesstyp == 1))));
        fprintf(ofp,"\tMean Capillary Flow = %.3f ± %.3f nl/mmin\n",mean(qq(find(vesstyp == 2))),stddev(qq(find(vesstyp == 2))));
        fprintf(ofp,"\tMean Venular Flow = %.3f ± %.3f nl/min\n",mean(qq(find(vesstyp == 3))),stddev(qq(find(vesstyp == 3))));
    }
    fprintf(ofp,"Mean Network Velocity = %.3f ± %.3f mm/s\n",mean(vel),stddev(vel));
    if (vesstyp.max() == 3) {
        fprintf(ofp,"\tMean Arteriolar Velocity = %.3f ± %.3f mm/s\n",mean(vel(find(vesstyp == 1))),stddev(vel(find(vesstyp == 1))));
        fprintf(ofp,"\tMean Capillary Velocity = %.3f ± %.3f mm/s\n",mean(vel(find(vesstyp == 2))),stddev(vel(find(vesstyp == 2))));
        fprintf(ofp,"\tMean Venular Velocity = %.3f ± %.3f mm/s\n",mean(vel(find(vesstyp == 3))),stddev(vel(find(vesstyp == 3))));
    }
    fprintf(ofp,"Wall Shear Stress = %.3f ± %.3f dyn/cm2\n",mean(tau),stddev(tau));
    if (vesstyp.max() == 3) {
        fprintf(ofp,"\tMean Arteriolar Wall Shear Stress = %.3f ± %.3f dyn/cm2\n",mean(tau(find(vesstyp == 1))),stddev(tau(find(vesstyp == 1))));
        fprintf(ofp,"\tMean Capillary Wall Shear Stress = %.3f ± %.3f dyn/cm2\n",mean(tau(find(vesstyp == 2))),stddev(tau(find(vesstyp == 2))));
        fprintf(ofp,"\tMean Venular Wall Shear Stress = %.3f ± %.3f dyn/cm2\n",mean(tau(find(vesstyp == 3))),stddev(tau(find(vesstyp == 3))));
    }
    fprintf(ofp,"Mean Vascular Segment Transit Time = %.3f ± %.3f s\n",mean(vstt),stddev(vstt));
    fprintf(ofp,"\nFry Parameters:\n");
    fprintf(ofp,"ktau = %f\n",ktau/2);
    fprintf(ofp,"kp = %f\n",kp);
    fprintf(ofp,"Target Pressure = %f mmHg\n",targPress);
    fprintf(ofp,"Target Wall Shear Stress = %f dyn/cm2\n",targStress);
    

    
    fclose(ofp);
    
}
