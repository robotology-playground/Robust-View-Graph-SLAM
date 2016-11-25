/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
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
