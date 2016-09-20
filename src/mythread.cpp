#include "mythread.h"
#include <iostream>

using namespace std;

bool MyThread::threadInit(){
    cout<<"Initializing thread number "<<id<<endl;
    return true;

}

void MyThread::run(){
    t->setFirstFrame(image1_cv);
    t->process(image2_cv);
}

void MyThread::threadRelease(){
    cout<<"Closing Thread number "<<id<<endl;
}


