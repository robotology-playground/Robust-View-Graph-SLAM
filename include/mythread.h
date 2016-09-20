#ifndef MYTHREAD_H
#define MYTHREAD_H

#include <yarp/os/Thread.h>
#include <Tracker.h>
#include <opencv2/features2d/features2d.hpp>
using namespace cv;

class MyThread : public yarp::os::Thread
{
    Tracker* t;
    int id;
    Mat image1_cv;
    Mat image2_cv;

public:
    MyThread(Ptr<Feature2D>& detector, Ptr<Feature2D>& descriptor, Ptr<DescriptorMatcher>& matcher,
             Mat& im1, Mat& im2,int num):yarp::os::Thread(){id=num;
                                                            t=new Tracker(detector,descriptor,matcher);
                                                            image1_cv=im1;
                                                            image2_cv=im2;
                                                           }
    ~MyThread(){
        delete t;
    }
protected:
    bool threadInit();
    void run();
    void threadRelease();

};

#endif // MYTHREAD_H
