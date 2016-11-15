/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threaddescriptor.h"


using namespace cv;

ThreadDescriptor::ThreadDescriptor(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut, cv::Ptr<cv::Feature2D> _descriptor)
    : vgSLAMThread(bufferIn, bufferOut),descriptor(_descriptor)
{

}

void ThreadDescriptor::run(){
    while(!interrupted) {
        SlamType data;
        yInfo()<<"ThreadDescriptor:Reading buffer";
        if(bufferIn->read(data)){
            data.descriptor= new Mat();
            descriptor->compute(*data.image, *data.feature, *data.descriptor);
            yInfo()<<"Descriptors:" << (int)data.descriptor->rows << ", " << (int)data.descriptor->cols;
            //bufferOut->write(data);
            data.free();
            countProcessed++;

        }
        else
            yInfo()<<"ThreadDescriptor has been interrupted on read";


    }
}
