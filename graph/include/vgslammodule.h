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
    yarp::os::BufferedPort<yarp::os::Bottle> encHeadPort, encTorsoPort;
    int imageL_start=0, imageR_start=0;
    bool first=true;
    ThreadFeature* threadFeatureR; //thread that reads from a buffer of R images and writes to a buffer of R features
    ThreadFeature* threadFeatureL; //thread that reads from a buffer of L images and writes to a buffer of L features
    ThreadDescriptor* threadDescriptorR; //thread that reads from a buffer of R features and writes to a prioritybuffer(*) of descriptors
    ThreadDescriptor* threadDescriptorL; //thread that reads from a buffer of L features and writes to a prioritybuffer(*) of descriptors
    ThreadMatching* threadMatching; //thread that reads from a prioritybuffer of descriptors, does matching and compute R&t

    vgSLAMBuffer<SlamType> bufferImageR;
    vgSLAMBuffer<SlamType> bufferImageL;
    vgSLAMBuffer<SlamType> bufferFeatureR;
    vgSLAMBuffer<SlamType> bufferFeatureL;
    vgSLAMPriorityBuffer<SlamType, SlamTypeComparison> bufferDescriptor;//prioritybuffer(*)
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
