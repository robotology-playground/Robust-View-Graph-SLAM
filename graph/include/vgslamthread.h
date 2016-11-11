/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email:nicolo.genesio@iit.it
 */

#ifndef VGSLAMTHREAD_H
#define VGSLAMTHREAD_H
#include "yarp/os/all.h"
#include "vgslambuffer.h"
template <class T1, class T2>
class vgSLAMThread : public yarp::os::Thread
{
    vgSLAMBuffer<T1> bufferIn;
    vgSLAMBuffer<T2> bufferOut;


public:
    vgSLAMThread() { }
    virtual void run() { }

};

#endif // VGSLAMTHREAD_H
