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


class ThreadMatching : public vgSLAMThread<SlamType,SlamType>
{
private:
    cv::Ptr<cv::DescriptorMatcher> matcher;
    SlamType *first;
    SlamType *second;
public:
    ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,cv::Ptr<cv::DescriptorMatcher> _matcher);
    virtual ~ThreadMatching();
    void run ();
protected:
    cv::Mat getProjMat(MatchesVector &matches, SlamType *data1, SlamType *data2);
    void getKinTransformationsToRoot(cv::Mat& ProjectionMatrix, SlamType *data);
    boost::dynamic_bitset<> remove_points_at_infinity(std::vector<cv::Point2f> matched1, std::vector<cv::Point2f> matched2, double thres);
};

#endif // THREADMATCHING_H
