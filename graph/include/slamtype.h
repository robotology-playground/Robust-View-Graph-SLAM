/*
 *  Created on: Nov 14, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */


#ifndef SLAMTYPE
#define SLAMTYPE

#include "opencvs.h"
#include "yarp/os/all.h"

typedef std::vector<cv::KeyPoint> KeyPointsVector;

class SlamType {
public:
    SlamType() : image(NULL), feature(NULL), descriptor(NULL),
    matching(NULL), stamp(NULL), relative(NULL) { }
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
    }
//    bool operator<(const SlamType& rhs)
//    {
//      return this->stamp->getTime() < rhs.stamp->getTime();
//    }

public:
    cv::Mat *image;
    KeyPointsVector *feature;
    cv::Mat *descriptor;
    cv::DMatch *matching;
    cv::Mat  *relative;
    yarp::os::Stamp *stamp;

};

/*
static bool  operator<(const SlamType& lhs, const SlamType& rhs)
{
    yDebug()<<__LINE__;
    return lhs.stamp->getCount() < rhs.stamp->getCount();
}
*/

class SlamTypeComparison
{
public:
    bool operator() (const SlamType& lhs, const SlamType& rhs) const
    {
      return lhs.stamp->getTime() < rhs.stamp->getTime();
    }
};

#endif // SLAMTYPE

