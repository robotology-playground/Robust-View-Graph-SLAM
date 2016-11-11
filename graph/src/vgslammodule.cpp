/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "vgslammodule.h"

vgSLAMModule::vgSLAMModule():nCams(-1)
{

}

vgSLAMModule::vgSLAMModule(int _nCams){
    nCams=_nCams;
}

bool vgSLAMModule::configure(yarp::os::ResourceFinder &rf){
    //Open ports
    yarp::os::Network yarp; // 64 bytes still reachable
    imageR_port.open("/vgSLAM/cam/left");
    imageL_port.open("/vgSLAM/cam/right");
    yarp.connect("/icub/cam/left", imageR_port.getName());
    yarp.connect("/icub/cam/right",imageL_port.getName());
    //start threads
    threadFeatureL.start();
    threadFeatureR.start();
    threadDescriptorL.start();
    threadDescriptorR.start();
    //threadMatching.start();

    return true;
}

bool vgSLAMModule::updateModule(){
    std::cout<<"Do something"<<std::endl;
}

bool vgSLAMModule::close(){
    imageR_port.close();
    imageL_port.close();
    threadFeatureL.stop();
    threadFeatureR.stop();
    threadDescriptorL.stop();
    threadDescriptorR.stop();
    //threadMatching.stop();
    return true;
}

bool vgSLAMModule::interruptModule(){
    this->close();
    imageR_port.interrupt();
    imageL_port.interrupt();
    return true;
}

