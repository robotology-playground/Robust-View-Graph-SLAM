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
    //We continue to read, process and write data while it is not interrupted and there is data to process.
    //The interrupted flag is setted to true in the function interrupt() of the parent class "vgSLAMThread".
    //interrupt() is called by thread.close() when we close the module.
    while(!interrupted) {
        SlamType data;
        yInfo()<<"ThreadDescriptor:Reading buffer";
        if(bufferIn->read(data)){
            data.descriptor= new Mat();
            //We compute descriptors and we write in the bufferOut, in this case the bufferOut is a priorityBuffer, it is shared
            //with the other ThreadDescriptor. It stores "data"(that contains information about images, features etc) ordered by
            //time.
            descriptor->compute(*data.image, *data.feature, *data.descriptor);
            yInfo()<<"Descriptors:" << (int)data.descriptor->rows << ", " << (int)data.descriptor->cols;
            bufferOut->write(data);

        }
        else
            yInfo()<<"ThreadDescriptor has been interrupted on read";


    }
}
