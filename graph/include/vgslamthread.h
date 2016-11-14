/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
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

    void interrupt() {
        interrupted = true;
        bufferIn->interrupt();
    }

    void close() {
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
