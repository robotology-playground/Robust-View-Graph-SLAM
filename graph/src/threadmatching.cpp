/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadmatching.h"

#include <stdio.h>

ThreadMatching::ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,
                               cv::Ptr<cv::DescriptorMatcher> _matcher): vgSLAMThread(bufferIn, bufferOut),matcher(_matcher){

}
void ThreadMatching::run (){
    while(!interrupted) {
        SlamType data;
        yInfo()<<"ThreadMatching:Reading buffer";
        if(bufferIn->read(data)){
            char str[255];
            //sprintf(str, "%.6f", data.stamp->getTime());
            //yDebug()<<"ThreadMatching:           "<<data.stamp->getCount();
            //bufferOut->write(data);
            data.free();
            countProcessed++;

        }
        else
            yInfo()<<"ThreadMatching has been interrupted on read";


    }

}
