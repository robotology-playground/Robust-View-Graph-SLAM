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
    matching(NULL), stamp(NULL) { }
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
        if(stamp){
            delete stamp;
            stamp=NULL;
        }
    }

public:
    cv::Mat *image;
    KeyPointsVector *feature;
    cv::Mat *descriptor;
    cv::DMatch *matching;
    yarp::os::Stamp *stamp;

};

#endif // SLAMTYPE

