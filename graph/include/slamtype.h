/**
 * @file slamtype.h
 * @brief This class all the types of pointers used in vgSLAM.
 * @detail It is also implemented the comparison for this kind of data type.
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#ifndef SLAMTYPE
#define SLAMTYPE

#include "opencvs.h"
#include "yarp/os/all.h"

typedef std::vector<cv::KeyPoint> KeyPointsVector;
typedef std::vector<cv::DMatch> MatchesVector;

//Class that includes all the types of data that threads will manage
class SlamType {
public:
    SlamType() : image(NULL), feature(NULL), descriptor(NULL),
    matching(NULL), stamp(NULL), relative(NULL),anglesHead(NULL), anglesTorso(NULL),
    right(true){ }
    virtual ~SlamType() {}
    void free() {
        if(image) {
            delete image;
            image = NULL;//prevent double call of free()
        }
        if(feature) {
            delete feature;
            feature = NULL;
        }
        if(descriptor) {
            delete descriptor;
            descriptor = NULL;
        }
        if(matching) {
            delete matching;
            matching = NULL;
        }
        if(relative) {
            delete relative;
            relative = NULL;
        }
        if(stamp){
            delete stamp;
            stamp = NULL;
        }
        if(anglesHead){
            delete anglesHead;
            anglesHead = NULL;
        }
        if(anglesTorso){
            delete anglesTorso;
            anglesTorso = NULL;
        }
    }
//    bool operator<(const SlamType& rhs)
//    {
//      return this->stamp->getTime() < rhs.stamp->getTime();
//    }

public:
    cv::Mat *image;
    std::vector<double> *anglesHead;
    std::vector<double> *anglesTorso;
    KeyPointsVector *feature;
    cv::Mat *descriptor;
    MatchesVector *matching;
    cv::Mat  *relative;
    bool right;
    yarp::os::Stamp *stamp;

};

/*
static bool  operator<(const SlamType& lhs, const SlamType& rhs)
{
    yDebug()<<__LINE__;
    return lhs.stamp->getCount() < rhs.stamp->getCount();
}
*/
// We define this class for using the priority_queue of the vgSLAMPriorityBuffer, in this way in this kind of buffer the elements
// are pushed ordered by the timestamp.
class SlamTypeComparison
{
public:
    bool operator() (const SlamType& lhs, const SlamType& rhs) const
    {
        //printf("\n%.6f\t%.6f\t%d\n\n", lhs.stamp->getTime(), rhs.stamp->getTime(), (lhs.stamp->getTime() < rhs.stamp->getTime()));
        return lhs.stamp->getTime() < rhs.stamp->getTime();
    }
};

#endif // SLAMTYPE

