#ifndef THREADMATCHING_H
#define THREADMATCHING_H

#include "vgslamthread.h"


class ThreadMatching : public vgSLAMThread<int,int>
{
public:
    ThreadMatching();
};

#endif // THREADMATCHING_H
