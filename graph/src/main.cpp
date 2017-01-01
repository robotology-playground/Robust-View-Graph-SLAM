/**
 * @file main.cpp
 * @brief A test file to demonstrate the multi-threading functionality of view-graph slam
 * @detail .
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */
#include "featureselector.h"
#include "yarp/os/all.h"
#include "yarp/sig/all.h"
#include "yarp/math/Math.h"
#include "Tracker.h"
#include "vgslammodule.h"
#include "stdio.h"
#include "iostream"

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;

int main (int argc, char** argv) {
    yarp::os::Network yarp;
    vgSLAMModule module(10);
    yarp::os::ResourceFinder rf;
    rf.setDefaultConfigFile("../../conf/vgSLAM.ini");
    rf.configure(argc,argv);
    //it calls .configure and then .updateModule
    module.runModule(rf);
    return 0;
}
