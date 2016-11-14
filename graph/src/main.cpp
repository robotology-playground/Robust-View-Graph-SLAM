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
    module.runModule(rf);

    return 0;

}
