/**
 * @file threaddescriptor.cpp
 * @brief Thread for computing descriptors.
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
