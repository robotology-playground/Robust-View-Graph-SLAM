/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadfeature.h"
#include <yarp/os/LogStream.h>
#include <yarp/os/Time.h>

using namespace yarp::os;

ThreadFeature::ThreadFeature(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut, cv::Ptr<cv::Feature2D> _detector)
    : vgSLAMThread(bufferIn, bufferOut),detector(_detector)
{

}

//ThreadFeature::ThreadFeature(cv::Ptr<cv::Feature2D> _detector):detector(_detector){
//}

void ThreadFeature::run(){
    while(!interrupted) {
        SlamType data;
        if(bufferIn->read(data))
            yInfo()<<"ThreadFeature read:";
        else
            yWarning()<<"ThreadFeature has been interrupted on read";
        //Time::delay(0.001);
        //data.free();//lo deve fare l'ultimo thread
        countProcessed++;
    }
}
