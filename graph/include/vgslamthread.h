/**
 * @file vglamthread.h
 * @brief Template class for defining the threads used in vgSLAM
 * @detail Each thread has two buffer, one as input and one as output.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @author Nicolo' Genesio
 * @email nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef VGSLAMTHREAD_H
#define VGSLAMTHREAD_H
#include "yarp/os/all.h"
#include "vgslambuffer.h"

template <class T1, class T2>
class vgSLAMThread : public yarp::os::Thread
{

public:
    vgSLAMThread(vgSLAMBuffer<T1> &bufferIn, vgSLAMBuffer<T2> &bufferOut)
    {
        interrupted = false;
        countProcessed = 0;
        vgSLAMThread::bufferIn = &bufferIn;
        vgSLAMThread::bufferOut = &bufferOut;
    }

    unsigned int getCountProcessed(){
        return countProcessed;
    }

    virtual void interrupt() {
        interrupted = true;
        bufferIn->interrupt();
    }

    virtual void close() {
        interrupt();
        Thread::stop();
        countProcessed = 0;
        interrupted = false;
    }

protected:
    unsigned int countProcessed;
    bool interrupted;
    vgSLAMBuffer<T1> *bufferIn;
    vgSLAMBuffer<T2> *bufferOut;

};

#endif // VGSLAMTHREAD_H
