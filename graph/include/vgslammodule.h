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

class vgSLAMModule : public yarp::os::RFModule
{
    int nCams;
    yarp::os::BufferedPort<yarp::sig::ImageOf<yarp::sig::PixelRgb> > imageR_port, imageL_port;
    ThreadFeature* threadFeatureR;
    ThreadFeature* threadFeatureL;
    ThreadDescriptor* threadDescriptorR;
    ThreadDescriptor* threadDescriptorL;
    ThreadMatching* threadMatching;
    vgSLAMBuffer<cv::Mat> bufferOut;

public:
    vgSLAMModule();
    vgSLAMModule(int _nCams);
   bool configure(yarp::os::ResourceFinder &rf);
   bool updateModule();
   bool interruptModule();
   bool close();
};

#endif // VGSLAMMODULE_H
