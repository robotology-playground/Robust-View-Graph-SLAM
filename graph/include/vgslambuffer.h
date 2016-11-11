/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email:nicolo.genesio@iit.it
 */

#ifndef VGSLAMBUFFER_H
#define VGSLAMBUFFER_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <yarp/os/all.h>
template <class T>
class vgSLAMBuffer
{
private:
    yarp::os::Mutex bufMutex;
    yarp::os::Event event;
    std::vector<T> array;
    bool interrupted;

public:
    vgSLAMBuffer() : interrupted(false) {}

    bool read(T& data) {
        bufMutex.lock();
        bool empty = array.empty();
        bufMutex.unlock();

        if(empty && !interrupted)
            event.wait();

        if(interrupted) {
            return false;
        }

        bufMutex.lock();
        data = array.front();
        bufMutex.unlock();
        return true;
    }

    bool write(T& data) {
        bufMutex.lock();
        array.push_back(data);
        bufMutex.unlock();
        event.signal();
        return true;
    }

    void interrupt() {
        bufMutex.lock();
        interrupted = true;
        bufMutex.unlock();
        event.signal();
    }
};

#endif // VGSLAMBUFFER_H
