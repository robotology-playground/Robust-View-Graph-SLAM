/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadfeature.h"

ThreadFeature::ThreadFeature()
{

}

ThreadFeature::ThreadFeature(cv::Ptr<cv::Feature2D> _detector):detector(_detector){

}
