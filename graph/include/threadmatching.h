/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef THREADMATCHING_H
#define THREADMATCHING_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"


class ThreadMatching : public vgSLAMThread<SlamType,SlamType>
{
private:
    cv::Ptr<cv::DescriptorMatcher> matcher;
public:
    ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::DescriptorMatcher> _matcher);
    void run ();
};

#endif // THREADMATCHING_H
