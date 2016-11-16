/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef VGSLAMMODULE_H
#define VGSLAMMODULE_H

#include "yarp/os/all.h"
#include "yarp/sig/all.h"
#include "vgslambuffer.h"
#include "threadfeature.h"
#include "threaddescriptor.h"
#include "threadmatching.h"
#include "opencv2/core.hpp"
#include "slamtype.h"

class vgSLAMModule : public yarp::os::RFModule
{
    bool configured;
    int nCams;
    yarp::os::BufferedPort<yarp::sig::ImageOf<yarp::sig::PixelRgb> > imageR_port, imageL_port;
    int imageL_start=0, imageR_start=0;
    bool first=true;
    ThreadFeature* threadFeatureR;
    ThreadFeature* threadFeatureL;
    ThreadDescriptor* threadDescriptorR;
    ThreadDescriptor* threadDescriptorL;
    ThreadMatching* threadMatching;

    vgSLAMBuffer<SlamType> bufferImageR;
    vgSLAMBuffer<SlamType> bufferImageL;
    vgSLAMBuffer<SlamType> bufferFeatureR;
    vgSLAMBuffer<SlamType> bufferFeatureL;
    vgSLAMPriorityBuffer<SlamType, SlamTypeComparison> bufferDescriptor;
    vgSLAMBuffer<SlamType> bufferMatching;


public:
    vgSLAMModule();
    vgSLAMModule(int _nCams);
   bool configure(yarp::os::ResourceFinder &rf);
   bool updateModule();
   double getPeriod();
   bool interruptModule();
   bool close();
};

#endif // VGSLAMMODULE_H
