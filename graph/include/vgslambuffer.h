/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef VGSLAMBUFFER_H
#define VGSLAMBUFFER_H

#include <stdio.h>
#include <iostream>
#include <queue>
#include <yarp/os/all.h>

template <class T>
class vgSLAMBuffer
{
private:
    yarp::os::Mutex bufMutex;
    yarp::os::Semaphore semArray;
    std::queue<T> array;
    bool interrupted;

public:
    vgSLAMBuffer() : interrupted(false), semArray(0){}

    bool read(T& data) {

        //wait
        semArray.wait();
        //yDebug()<<"read() : unwait";

        if(interrupted) {
            return false;
        }

        // read
        bufMutex.lock();
        data = array.front();
        array.pop();
        bufMutex.unlock();
        return true;
    }

    bool write(T& data) {
        //write
        bufMutex.lock();
        array.push(data);
        //yDebug()<<"write() : signal";
        bufMutex.unlock();
        //signal
        semArray.post();
        return true;
    }

    void interrupt() {
        bufMutex.lock();
        interrupted = true;
        bufMutex.unlock();
        semArray.post();
    }
    void clear(){
        bufMutex.lock();
        array.clear();
        bufMutex.unlock();
    }
};

#endif // VGSLAMBUFFER_H
