/**
 * @file threaddescriptor.h
 * @brief Thread for computing descriptors.
 * @detail It inerhits from the template vgSLAMThread
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef THREADDESCRIPTOR_H
#define THREADDESCRIPTOR_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"

class ThreadDescriptor : public vgSLAMThread<SlamType, SlamType>
{
private:
    cv::Ptr<cv::Feature2D> descriptor;
public:
    ThreadDescriptor(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::Feature2D> _descriptor);
    void run ();
};

#endif // THREADDESCRIPTOR_H
