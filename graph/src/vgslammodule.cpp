/**
 * @file vgslammodule.h
 * @brief Main core of vgSLAM.
 * @detail It defines, configures, starts and stops all the thread of vgSLAM.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */
#include <yarp/os/LogStream.h>
#include "vgslammodule.h"
#include "featureselector.h"


using namespace yarp::os;
using namespace yarp::sig;
using namespace cv;

vgSLAMModule::vgSLAMModule():nCams(-1),  bufferDescriptorR(),bufferDescriptorL()
{

}

vgSLAMModule::vgSLAMModule(int _nCams) :  bufferDescriptorR(),bufferDescriptorL()
{
    nCams=_nCams;
    configured = false;
}

//double anglesHead[6], anglesTorso[3];
//Bottle *headBottle, *torsoBottle;
//BufferedPort<Bottle> encHeadPort, encTorsoPort;
bool vgSLAMModule::configure(ResourceFinder &rf){
    //Open ports
    imageR_port.open("/vgSLAM/cam/left");
    imageL_port.open("/vgSLAM/cam/right");
    encHeadPort.open("/vgSLAM/head/state:i");
    encTorsoPort.open("/vgSLAM/torso/state:i");
    //Connect ports
    bool ret = NetworkBase::connect("/icub/cam/left", imageR_port.getName());
    ret &= NetworkBase::connect("/icub/cam/right",imageL_port.getName());
    ret &= NetworkBase::connect("/icub/head/state:o",encHeadPort.getName());
    ret &= NetworkBase::connect("/icub/torso/state:o",encTorsoPort.getName());
    if(!ret) {
        yError()<<"Could not connect to some of the ports";
        return configured;
    }
    //configure threads and select detector, descriptor and matcher
    Ptr<Feature2D> detector;
    Ptr<Feature2D> descriptor;
    Ptr<DescriptorMatcher> matcher;
    FeatureSelector selector(rf);
    selector.process(detector,descriptor,matcher);
    threadFeatureL = new ThreadFeature(bufferImageL, bufferFeatureL,detector);
    threadFeatureR = new ThreadFeature(bufferImageR, bufferFeatureR,detector);
    threadDescriptorL = new ThreadDescriptor(bufferFeatureL,bufferDescriptorL,descriptor);
    threadDescriptorR = new ThreadDescriptor(bufferFeatureR,bufferDescriptorR,descriptor);
    threadMatching=new ThreadMatching(bufferDescriptorL,bufferDescriptorR,bufferMatching,matcher);

    //start threads
    threadFeatureL->start();
    threadFeatureR->start();
    threadDescriptorL->start();
    threadDescriptorR->start();
    threadMatching->start();

    configured = true;
    return configured;
}

//It is called while it returns true
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
    Bottle  *headBottle=encHeadPort.read();
    Bottle  *torsoBottle=encTorsoPort.read();
    dataL.anglesHead=new std::vector<double>(6);
    dataR.anglesHead=new std::vector<double>(6);
    dataL.anglesTorso=new std::vector<double>(6);
    dataR.anglesTorso=new std::vector<double>(6);

    dataL.anglesHead->at(0)=headBottle->get(0).asDouble();
    dataL.anglesHead->at(1)=headBottle->get(1).asDouble();
    dataL.anglesHead->at(2)=headBottle->get(2).asDouble();
    dataL.anglesHead->at(3)=headBottle->get(3).asDouble();
    dataL.anglesHead->at(4)=headBottle->get(4).asDouble();
    dataL.anglesHead->at(5)=headBottle->get(5).asDouble();

    dataR.anglesHead->at(0)=headBottle->get(0).asDouble();
    dataR.anglesHead->at(1)=headBottle->get(1).asDouble();
    dataR.anglesHead->at(2)=headBottle->get(2).asDouble();
    dataR.anglesHead->at(3)=headBottle->get(3).asDouble();
    dataR.anglesHead->at(4)=headBottle->get(4).asDouble();
    dataR.anglesHead->at(5)=headBottle->get(5).asDouble();

    dataL.anglesTorso->at(0)=torsoBottle->get(0).asDouble();
    dataL.anglesTorso->at(1)=torsoBottle->get(1).asDouble();
    dataL.anglesTorso->at(2)=torsoBottle->get(2).asDouble();
    dataL.anglesTorso->at(3)=torsoBottle->get(3).asDouble();
    dataL.anglesTorso->at(4)=torsoBottle->get(4).asDouble();
    dataL.anglesTorso->at(5)=torsoBottle->get(5).asDouble();

    dataR.anglesTorso->at(0)=torsoBottle->get(0).asDouble();
    dataR.anglesTorso->at(1)=torsoBottle->get(1).asDouble();
    dataR.anglesTorso->at(2)=torsoBottle->get(2).asDouble();
    dataR.anglesTorso->at(3)=torsoBottle->get(3).asDouble();
    dataR.anglesTorso->at(4)=torsoBottle->get(4).asDouble();
    dataR.anglesTorso->at(5)=torsoBottle->get(5).asDouble();
// Tested
//    yDebug()<<"TorsoR:"<<dataR.anglesTorso->at(0)<<dataR.anglesTorso->at(1)<<dataR.anglesTorso->at(2)<<dataR.anglesTorso->at(3)<<dataR.anglesTorso->at(4)<<dataR.anglesTorso->at(5);
//    yDebug()<<"TorsoL:"<<dataL.anglesTorso->at(0)<<dataL.anglesTorso->at(1)<<dataL.anglesTorso->at(2)<<dataL.anglesTorso->at(3)<<dataL.anglesTorso->at(4)<<dataL.anglesTorso->at(5);
//    yDebug()<<"HeadR:"<<dataR.anglesHead->at(0)<<dataR.anglesHead->at(1)<<dataR.anglesHead->at(2)<<dataR.anglesHead->at(3)<<dataR.anglesHead->at(4)<<dataR.anglesHead->at(5);
//    yDebug()<<"HeadL:"<<dataL.anglesHead->at(0)<<dataL.anglesHead->at(1)<<dataL.anglesHead->at(2)<<dataL.anglesHead->at(3)<<dataL.anglesHead->at(4)<<dataL.anglesTorso->at(5);

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
// syncronization: ask to Alberto for syncronizer
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
        dataL.right=false;
        dataR.image = new Mat(imageR_cv);
        dataR.right=true;

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
//Once all the images are acquired the module will be closed, but first it waits that all the threads
//finish their work. We are checking how many frames has been processed by the last thread(matching)
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
    //This method is defined in order to call updateModule in an infinite loop while it returns true.
    return 0.0;
}

bool vgSLAMModule::interruptModule(){
    imageR_port.interrupt();
    imageL_port.interrupt();
    return true;
}

