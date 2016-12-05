/**
 * @file threadfeatures.cpp
 * @brief Thread for features detection.
 * @detail It inerhits from the template vgSLAMThread
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @author Nicolo' Genesio
 * @email nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
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
