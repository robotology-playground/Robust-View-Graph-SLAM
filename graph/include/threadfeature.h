/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef THREADFEATURE_H
#define THREADFEATURE_H

#include "vgslamthread.h"
#include "opencvs.h"

class ThreadFeature : public vgSLAMThread<int,int>
{
private:
    cv::Ptr<cv::Feature2D> detector;
public:
    ThreadFeature();
    ThreadFeature(cv::Ptr<cv::Feature2D> _detector);
};

#endif // THREADFEATURE_H
