/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadfeature.h"
#include <yarp/os/LogStream.h>
#include <yarp/os/Time.h>


using namespace yarp::os;
using namespace cv;

ThreadFeature::ThreadFeature(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut, cv::Ptr<cv::Feature2D> _detector)
    : vgSLAMThread(bufferIn, bufferOut),detector(_detector)
{

}


void ThreadFeature::run(){
    while(!interrupted) {
        SlamType data;
        yInfo()<<"ThreadFeature:Reading buffer";
        if(bufferIn->read(data)){
            data.feature = new KeyPointsVector;
            detector->detect(*data.image, *data.feature, noArray());
            yInfo()<<"ThreadFeature read: Number of keypoints="<<data.feature->size();
            bufferOut->write(data);
        }
        else
            yInfo()<<"ThreadFeature has been interrupted on read";

    }
}
