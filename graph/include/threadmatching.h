/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email:nicolo.genesio@iit.it
 */

#ifndef THREADMATCHING_H
#define THREADMATCHING_H

#include "vgslamthread.h"


class ThreadMatching : public vgSLAMThread<int,int>
{
public:
    ThreadMatching();
};

#endif // THREADMATCHING_H
