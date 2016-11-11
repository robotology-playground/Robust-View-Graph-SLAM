/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email:nicolo.genesio@iit.it
 */

#ifndef THREADDESCRIPTOR_H
#define THREADDESCRIPTOR_H

#include "vgslamthread.h"


class ThreadDescriptor : public vgSLAMThread<int,int>
{
public:
    ThreadDescriptor();
};

#endif // THREADDESCRIPTOR_H
