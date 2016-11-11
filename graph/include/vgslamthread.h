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
    vgSLAMBuffer<T1> bufferIn;
    vgSLAMBuffer<T2> bufferOut;


public:
    vgSLAMThread() { }
    virtual void run() {
        T1 dataIn;
        T2 dataOut;
        bufferIn.read(dataIn);
        dataOut=dataOut+dataIn;
        bufferOut.write(dataOut);
    }
    virtual void stop() {
        bufferIn.clear();
        bufferOut.clear();
    }
//    virtual void run()=0;
//    virtual void process()=0; dovranno essere cosi', almeno uno no puo' usare la classe senza implementarla.
//    virtual void stop()=0;


};

#endif // VGSLAMTHREAD_H
