#ifndef VGSLAMMODULE_H
#define VGSLAMMODULE_H

#include "yarp/os/all.h"
#include "vgslambuffer.h"
#include "threadfeature.h"
#include "threaddescriptor.h"
#include "threadmatching.h"

class vgSLAMModule : public yarp::os::RFModule
{
    ThreadFeature threadFeatureR;
    ThreadFeature threadFeatureL;
    ThreadDescriptor threadDescriptorR;
    ThreadDescriptor threadDescriptorL;
    ThreadMatching threadMatching;
    vgSLAMBuffer<int> bufferOut;

public:
    vgSLAMModule();
};

#endif // VGSLAMMODULE_H
