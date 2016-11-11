/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#ifndef THREADMATCHING_H
#define THREADMATCHING_H

#include "vgslamthread.h"
#include "opencvs.h"

class ThreadMatching : public vgSLAMThread<int,int>
{
private:
    cv::Ptr<cv::DescriptorMatcher> matcher;
public:
    ThreadMatching();
    ThreadMatching(cv::Ptr<cv::DescriptorMatcher> _matcher);
};

#endif // THREADMATCHING_H
