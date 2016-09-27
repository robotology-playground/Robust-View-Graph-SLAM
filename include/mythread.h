#ifndef MYTHREAD_H
#define MYTHREAD_H

#include <yarp/os/Thread.h>
#include <Tracker.h>
#include <opencv2/features2d/features2d.hpp>
using namespace cv;

class MyThread : public yarp::os::Thread
{
    Tracker t;
    int id;
    Mat image_cv;
    Mat* ProjMat;

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
