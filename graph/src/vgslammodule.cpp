/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "vgslammodule.h"
#include "featureselector.h"

vgSLAMModule::vgSLAMModule():nCams(-1)
{

}

vgSLAMModule::vgSLAMModule(int _nCams){
    nCams=_nCams;
}

bool vgSLAMModule::configure(yarp::os::ResourceFinder &rf){
    //Open ports
    imageR_port.open("/vgSLAM/cam/left");
    imageL_port.open("/vgSLAM/cam/right");
    yarp::os::NetworkBase::connect("/icub/cam/left", imageR_port.getName());
    yarp::os::NetworkBase::connect("/icub/cam/right",imageL_port.getName());

    //configure threads
    cv::Ptr<cv::Feature2D> detector;
    cv::Ptr<cv::Feature2D> descriptor;
    cv::Ptr<cv::DescriptorMatcher> matcher;
    FeatureSelector selector(rf);
    selector.process(detector,descriptor,matcher);
    threadFeatureL=new ThreadFeature(detector);
    threadFeatureR=new ThreadFeature(detector);
//    threadDescriptorL=new ThreadDescriptor(descriptor);
//    threadDescriptorR=new ThreadDescriptor(descriptor);
//    threadMatching=new ThreadMatching(matcher);

    //start threads
    threadFeatureL->start();
    threadFeatureR->start();
//    threadDescriptorL->start();
//    threadDescriptorR->start();
    //threadMatching->start();

    return true;
}

bool vgSLAMModule::updateModule(){
    std::cout<<"Do something"<<std::endl;
}

bool vgSLAMModule::close(){
    //stop threads
    threadFeatureL->stop();
    threadFeatureR->stop();
//    threadDescriptorL->stop();
//    threadDescriptorR->stop();
//    threadMatching->stop();

    //close ports
    imageR_port.close();
    imageL_port.close();

    //deallocate memory
    delete threadFeatureL;
    delete threadFeatureR;
//    delete threadDescriptorL;
//    delete threadDescriptorR;
//    delete threadMatching;
    return true;
}

bool vgSLAMModule::interruptModule(){
    this->close();
    imageR_port.interrupt();
    imageL_port.interrupt();
    return true;
}

