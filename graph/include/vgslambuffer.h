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
#include <functional>

#include <slamtype.h>

template <class T>
class vgSLAMBuffer
{
protected:
    yarp::os::Mutex bufMutex;
    yarp::os::Semaphore semArray;
    std::queue<T> array;
    bool interrupted;

public:
    vgSLAMBuffer() : interrupted(false), semArray(0) {}

    virtual bool read(T& data) {

        //wait
        semArray.wait();
        //yDebug()<<"read() : unwait";

        if(interrupted) {
            return false;
        }

        // read
        bufMutex.lock();
        data=array.front();
        array.pop();
        bufMutex.unlock();
        return true;
    }

    virtual bool write(T& data) {
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

};



template <class T, class TComparison>
class vgSLAMPriorityBuffer : public vgSLAMBuffer<T> {

    typedef std::priority_queue<T, std::vector<T>, TComparison> mypq_type;

private:
    mypq_type p_array;

public:
    vgSLAMPriorityBuffer() : p_array(TComparison()) {}

    virtual bool read(T& data) {

        //wait
        vgSLAMBuffer<T>::semArray.wait();
        if(vgSLAMBuffer<T>::interrupted) {
            return false;
        }

        // read
        vgSLAMBuffer<T>::bufMutex.lock();
        data = p_array.top();
//        yDebug()<<"Priority access";//ok
        p_array.pop();
        vgSLAMBuffer<T>::bufMutex.unlock();
        return true;
    }

    virtual bool write(T& data) {
        //write
        vgSLAMBuffer<T>::bufMutex.lock();
        p_array.push(data);
        vgSLAMBuffer<T>::bufMutex.unlock();
        //signal
        vgSLAMBuffer<T>::semArray.post();
        return true;
    }
};

#endif // VGSLAMBUFFER_H
