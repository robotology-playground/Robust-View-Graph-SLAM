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
    //We continue to read, process and write data while it is not interrupted and there is data to process.
    //The interrupted flag is setted to true in the function interrupt() of the parent class "vgSLAMThread".
    //interrupt() is called by thread.close() when we close the module.
    while(!interrupted) {
        SlamType data;
        yInfo()<<"ThreadFeature:Reading buffer";
        if(bufferIn->read(data)){
            data.feature = new KeyPointsVector;
            //We compute features and we write in the bufferOut
            detector->detect(*data.image, *data.feature, noArray());
            yInfo()<<"ThreadFeature read: Number of keypoints="<<data.feature->size();
            bufferOut->write(data);
        }
        else
            yInfo()<<"ThreadFeature has been interrupted on read";

    }
}
