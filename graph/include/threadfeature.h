#ifndef THREADFEATURE_H
#define THREADFEATURE_H

#include "vgslamthread.h"


class ThreadFeature : public vgSLAMThread<int,int>
{
public:
    ThreadFeature();
};

#endif // THREADFEATURE_H
