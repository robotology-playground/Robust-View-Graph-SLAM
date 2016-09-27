/*
 * mythread.cpp
 *
 *  Created on: Sept 21, 2016
 *      Author: Nicolo' Genesio
 * This class inherits from yarp::os::Thread class.
 * It has been designed to compute feature and do matching with the reference frame in a
 * multithreading context.
 * The method ThreadInit(), run() and threadRelease() has been overridden from the Parent class.
 */

#include "mythread.h"
#include <iostream>
#include <yarp/os/LogStream.h>

using namespace std;

bool MyThread::threadInit(){
    cout<<"Initializing thread number "<<id<<endl;
    return true;

}

void MyThread::run(){
    //yInfo()<<"Start run";

    t.process(image_cv,ProjMat);
    //yInfo()<<"End run";
}

void MyThread::threadRelease(){
    cout<<"Closing Thread number "<<id<<endl;
}
