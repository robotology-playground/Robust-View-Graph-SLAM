/*
 *  Created on: Nov 14, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */


#ifndef SLAMTYPE
#define SLAMTYPE

#include "opencvs.h"

class SlamType {
public:
    SlamType() : image(NULL), feature(NULL) { }
    virtual ~SlamType() { free(); }
    void free() {
        if(image)
            delete image;
        if(feature)
            delete feature;
    }

public:
    cv::Mat *image;
    cv::Point2f *feature;
};

#endif // SLAMTYPE

