/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadmatching.h"

ThreadMatching::ThreadMatching()
{

}

ThreadMatching::ThreadMatching(cv::Ptr<cv::DescriptorMatcher> _matcher):matcher(_matcher){

}
