/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include <yarp/os/LogStream.h>
#include "vgslammodule.h"
#include "featureselector.h"


using namespace yarp::os;
using namespace yarp::sig;
using namespace  cv;

vgSLAMModule::vgSLAMModule():nCams(-1),  bufferDescriptor()
{

}

vgSLAMModule::vgSLAMModule(int _nCams) :  bufferDescriptor(){
    nCams=_nCams;
    configured = false;
}

bool vgSLAMModule::configure(ResourceFinder &rf){
    //Open ports
    imageR_port.open("/vgSLAM/cam/left");
    imageL_port.open("/vgSLAM/cam/right");
    bool ret = NetworkBase::connect("/icub/cam/left", imageR_port.getName());
    ret &= NetworkBase::connect("/icub/cam/right",imageL_port.getName());
    if(!ret) {
        yError()<<"Could not connect to some of the ports";
        return configured;
    }
    //configure threads
    Ptr<Feature2D> detector;
    Ptr<Feature2D> descriptor;
    Ptr<DescriptorMatcher> matcher;
    FeatureSelector selector(rf);
    selector.process(detector,descriptor,matcher);
    threadFeatureL = new ThreadFeature(bufferImageL, bufferFeatureL,detector);
    threadFeatureR = new ThreadFeature(bufferImageR, bufferFeatureR,detector);
    threadDescriptorL = new ThreadDescriptor(bufferFeatureL,bufferDescriptor,descriptor);
    threadDescriptorR = new ThreadDescriptor(bufferFeatureR,bufferDescriptor,descriptor);
    threadMatching=new ThreadMatching(bufferDescriptor,bufferMatching,matcher);

    //start threads
    threadFeatureL->start();
    threadFeatureR->start();
    threadDescriptorL->start();
    threadDescriptorR->start();
    threadMatching->start();

    configured = true;
    return configured;
}

bool vgSLAMModule::updateModule(){
    static int count = nCams;
    if(count<=0) {
        yDebug()<<"finishing acquisition...";
        return false;
    }

    // read port
    SlamType dataL,dataR;

    yInfo() <<"vgSLAMModule:acquiring images...";
    ImageOf<PixelRgb> *imageL_yarp = imageL_port.read();
    ImageOf<PixelRgb> *imageR_yarp = imageR_port.read();

    if (!imageL_yarp || !imageR_yarp) {
        configured = false;
        return false;
    }

    yInfo() <<"vgSLAMModule:acquiring images [done]";

    Stamp sL,sR;
    if(imageL_port.getEnvelope(sL) && imageR_port.getEnvelope(sR)){
        if(first){
            imageL_start = sL.getCount();
            imageR_start = sR.getCount();
            first=false;
        }
// syncronization ask to Alberto for syncronizer
//        if(abs(sL.getCount()-imageL_start)>2 || abs(sR.getCount()-imageR_start)>2
//                || fabs((sL.getTime())-(sR.getTime()))>0.06){//0.03 is the half deltat (30hz) but it is more likely 15 hz
//            imageL_start = sL.getCount();
//            imageR_start = sR.getCount();
//            yWarning()<<"Left-Right de-synchronized, time difference:"<<fabs((sL.getTime())-(sR.getTime()));
//            return true;
//        }

        /* the images */
        Mat imageL_cv = cvarrToMat(static_cast<IplImage*>(imageL_yarp->getIplImage()));
        cvtColor(imageL_cv, imageL_cv, CV_RGB2BGR);
        cvtColor(imageL_cv, imageL_cv, COLOR_BGR2GRAY);
        Mat imageR_cv = cvarrToMat(static_cast<IplImage*>(imageR_yarp->getIplImage()));
        cvtColor(imageR_cv, imageR_cv, CV_RGB2BGR);
        cvtColor(imageR_cv, imageR_cv, COLOR_BGR2GRAY);

        dataL.image = new Mat(imageL_cv);
        dataR.image = new Mat(imageR_cv);

        dataL.stamp =new Stamp(sL);
        dataR.stamp =new Stamp(sR);

        //imshow("prova",dataL.image);
        count -= 2;

        yInfo() <<"vgSLAMModule:writing images to buffers...";
        bufferImageL.write(dataL);
        bufferImageR.write(dataR);

    }
    return true;
}

bool vgSLAMModule::close(){

    if(configured) {
        yInfo()<<"Waiting for worker thread...";
        while(threadMatching->getCountProcessed() < (nCams-3))//3->thickness
            yarp::os::Time::delay(0.5);
    }
    //stop threads
    threadFeatureL->close();
    threadFeatureR->close();
    threadDescriptorL->close();
    threadDescriptorR->close();
    threadMatching->close();

    //close ports
    imageR_port.close();
    imageL_port.close();

    //deallocate memory
    delete threadFeatureL;
    delete threadFeatureR;
    delete threadDescriptorL;
    delete threadDescriptorR;
    delete threadMatching;
    return true;
}

double vgSLAMModule::getPeriod(){
    return 0.0;//fa un loop infinito
}

bool vgSLAMModule::interruptModule(){
    imageR_port.interrupt();
    imageL_port.interrupt();
    return true;
}

