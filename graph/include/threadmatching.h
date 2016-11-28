/**
 * @file threadmatching.h
 * @brief Thread for feature matching.
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
#ifndef THREADMATCHING_H
#define THREADMATCHING_H

#include "vgslamthread.h"
#include "opencvs.h"
#include "slamtype.h"
#include "boost/dynamic_bitset.hpp"
#include "Sampling.h"

class ThreadMatching : public vgSLAMThread<SlamType,SlamType>
{

private:
    vgSLAMBuffer<SlamType> *bufferIn1;
    vgSLAMBuffer<SlamType> *bufferIn2;
    cv::Ptr<cv::DescriptorMatcher> matcher;
    SlamType *first;
    SlamType *second;
    Sampling sampler;//for testing
public:
    ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn1,vgSLAMBuffer<SlamType> &bufferIn2, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::DescriptorMatcher> _matcher);
    virtual ~ThreadMatching();
    void interrupt();
    void run ();
protected:
    cv::Mat getProjMat(MatchesVector &matches, SlamType *data1, SlamType *data2);
    void getKinTransformationsToRoot(cv::Mat& ProjectionMatrix, SlamType *data);
    boost::dynamic_bitset<> remove_points_at_infinity(std::vector<cv::Point2f> matched1, std::vector<cv::Point2f> matched2, double thres);
};

#endif // THREADMATCHING_H
