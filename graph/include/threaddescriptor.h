/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef THREADDESCRIPTOR_H
#define THREADDESCRIPTOR_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"

class ThreadDescriptor : public vgSLAMThread<SlamType, SlamType>
{
private:
    cv::Ptr<cv::Feature2D> descriptor;
public:
    ThreadDescriptor(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::Feature2D> _descriptor);
    void run ();
};

#endif // THREADDESCRIPTOR_H
