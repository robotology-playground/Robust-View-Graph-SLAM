/*
 * mythread.h
 *
 *  Created on: Sept 21, 2016
 *      Author: Nicolo' Genesio
 * This class inherits from yarp::os::Thread class.
 * It has been designed to compute feature and do matching with the reference frame in a
 * multithreading context.
 * The method ThreadInit(), run() and threadRelease() has been overridden from the Parent class.
 */
#ifndef MYTHREAD_H
#define MYTHREAD_H

#include <yarp/os/Thread.h>
#include <Tracker.h>
#include <opencv2/features2d/features2d.hpp>
using namespace cv;

class MyThread : public yarp::os::Thread
{
    Tracker t;//Tracker used to compute and matche features
    int id;//number of the thread
    Mat image_cv;//image on which compute features
    Mat* ProjMat;//output of the processing, a projection matrix 3x4

public:
    MyThread(Tracker tracker, Mat& im, Mat& _ProjMat,int num):yarp::os::Thread(){id=num;
                                                            t=tracker;
                                                            image_cv=im;
                                                            ProjMat=&_ProjMat;
                                                           }
protected:
    bool threadInit();
    void run();
    void threadRelease();

};

#endif // MYTHREAD_H
