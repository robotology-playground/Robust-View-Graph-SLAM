/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email:nicolo.genesio@iit.it
 */

#ifndef THREADFEATURE_H
#define THREADFEATURE_H

#include "vgslamthread.h"


class ThreadFeature : public vgSLAMThread<int,int>
{
public:
    ThreadFeature();
};

#endif // THREADFEATURE_H
