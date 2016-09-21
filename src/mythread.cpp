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
    t.process(image_cv);
    //yInfo()<<"End run";
}

void MyThread::threadRelease(){
    cout<<"Closing Thread number "<<id<<endl;
}
