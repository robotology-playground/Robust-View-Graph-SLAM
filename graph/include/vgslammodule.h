/**
 * @file vgslammodule.h
 * @brief Main core of vgSLAM.
 * @detail It defines, configures, starts and stops all the thread of vgSLAM.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @author Nicolo' Genesio
 * @email nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
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
    vgSLAMBuffer<SlamType> bufferDescriptorR;
    vgSLAMBuffer<SlamType> bufferDescriptorL;
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
