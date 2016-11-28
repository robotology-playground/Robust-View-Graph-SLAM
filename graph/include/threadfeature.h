/**
 * @file threadfeatures.h
 * @brief Thread for features detection.
 * @detail It inerhits from the template vgSLAMThread
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @author Nicolo' Genesio
 * @email nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef THREADFEATURE_H
#define THREADFEATURE_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"

class ThreadFeature : public vgSLAMThread<SlamType, SlamType>
{
private:
    cv::Ptr<cv::Feature2D> detector;

public:
    ThreadFeature(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::Feature2D> _detector);
    void run ();
};

#endif // THREADFEATURE_H
