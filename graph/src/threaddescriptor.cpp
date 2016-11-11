/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threaddescriptor.h"

ThreadDescriptor::ThreadDescriptor()
{

}
ThreadDescriptor::ThreadDescriptor(cv::Ptr<cv::Feature2D> _descriptor):descriptor(_descriptor){

}

