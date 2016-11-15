/*
 *  Created on: Nov 14, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */


#ifndef SLAMTYPE
#define SLAMTYPE

#include "opencvs.h"

typedef std::vector<cv::KeyPoint> KeyPointsVector;

class SlamType {
public:
    SlamType() : image(NULL), feature(NULL), descriptor(NULL),
    matching(NULL) { }
    virtual ~SlamType() {}
    void free() {
        if(image) {
            delete image;
            image = NULL;
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
    }

public:
    cv::Mat *image;
    KeyPointsVector *feature;
    cv::Mat *descriptor;
    cv::DMatch *matching;
};

#endif // SLAMTYPE

