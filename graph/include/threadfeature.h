/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef THREADFEATURE_H
#define THREADFEATURE_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"

class ThreadFeature : public vgSLAMThread<SlamType, SlamType>
{
private:
    cv::Ptr<cv::Feature2D> detector;

public:
    ThreadFeature(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::Feature2D> _detector);

    //ThreadFeature(cv::Ptr<cv::Feature2D> _detector);
    void run ();
};

#endif // THREADFEATURE_H
