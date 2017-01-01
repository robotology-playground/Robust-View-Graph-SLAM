/**
 * @file vgslambuffer.h
 * @brief  Buffers that can be handled from multiple thread avoiding race condition.
 * @detail There are two kind of buffers, one that contains a queue as collection and the other
 * that use the std::priority_queue for storing in function of the time.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
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
    yarp::os::Mutex bufMutex;//mutex to lock the array, for avoiding race condition
    yarp::os::Semaphore semArray;
    std::queue<T> array;
    bool interrupted;

public:
    vgSLAMBuffer() : interrupted(false), semArray(0) {}

    virtual bool read(T& data) {

        //wait for new data
        /*
         * Decrement the counter, even if we must wait to do that.  If the counter
         * would decrement below zero, the calling thread must stop and
         * wait for another thread to call Semaphore::post on this semaphore.
         */
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
        //signal that a new data is ready
        /*
         * Increment the counter.  If another thread is waiting to decrement the
         * counter, it is woken up.
         */
        semArray.post();
        return true;
    }

    void interrupt() {
        //it is called in vgSLAMThread::interrupt() that is called when we close the module and then the
        //threads.
        bufMutex.lock();
        interrupted = true;
        bufMutex.unlock();
        semArray.post();
    }

};



template <class T, class TComparison>
class vgSLAMPriorityBuffer : public vgSLAMBuffer<T> {
    //the priority_queue needs the type, and the operator() to do the comparison necessary for ordering the queue.
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
