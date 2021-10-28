//
//  Header.h
//  Vascular Flows 2.0
//
//  Created by Paul Sweeney on 05/06/2015.
//  Copyright (c) 2015 Paul Sweeney - University College London. All rights reserved.
//

#ifndef Vascular_Flows_2_0_Header_h
#define Vascular_Flows_2_0_Header_h

extern void input();
extern void setup_arrays();
extern void analyse_net();
extern void rheolParams();
extern double gasdev(long *idum);
extern double ran1(long *idum);
extern void check_dir(const int &i);

#endif
