/*
 * featureselector.h
 *
 *  Created on: Sept 7, 2016
 *      Author: Nicolo' Genesio
 * This class is designed to interface with the user through configuration files(.ini).
 * The user can select the method for the detector,descriptor and matcher in conf/vgSLAM.ini
 * and then tune the parameters of the method in the respective ini files contained in the conf
 * folder.
 */

#ifndef FEATURESELECTOR_H
#define FEATURESELECTOR_H

#define    vgSLAM_KAZE    0 //det
#define    vgSLAM_FAST    1 //det
#define    vgSLAM_SIFT    2 //det & desc
#define    vgSLAM_GFTT    3 //det
#define    vgSLAM_SURF    4 //det & desc 32f
#define    vgSLAM_BRIEF   5 //desc 8U
#define    vgSLAM_ORB     6 //det & desc 8U
#define    vgSLAM_BRISK   7 //det & desc 8U
#define    vgSLAM_FREAK   8 //desc 8U
#define    vgSLAM_FLANN   9 //matcher  works with 32f desc
#define    vgSLAM_BRUTEFORCEL2   10 //matcher
#define    vgSLAM_BRUTEFORCEL1 11 //matcher
#define    vgSLAM_BRUTEFORCEHAMMING   12//matcher works with 8U desc
#define    vgSLAM_AKAZE    13 //det

#include "yarp/os/all.h"
#include "opencvs.h"
#include "cpps.h"


class FeatureSelector
{
    yarp::os::ResourceFinder rf;


public:
    FeatureSelector();
    FeatureSelector(yarp::os::ResourceFinder _rf);
    bool process(cv::Ptr<cv::Feature2D> &detector, cv::Ptr<cv::Feature2D> &descriptor, cv::Ptr<cv::DescriptorMatcher> &matcher);
protected:
    int  parseMap(std::string value, std::map<std::string, int> &m);
    bool checkDetector(std::string str);//check if the method indicated in vgSLAM.ini is one of the available ones.
    bool checkDescriptor(std::string str);//check if the method indicated in vgSLAM.ini is one of the available ones.
    bool checkMatcher(std::string str, int _descID);//check if the method indicated in vgSLAM.ini is one of the available ones.
    void switcher(int detID, int descID, int matchID,
                  cv::Ptr<cv::Feature2D> &detector, cv::Ptr<cv::Feature2D> &descriptor,cv::Ptr<cv::DescriptorMatcher> &matcher);
    //given the id for the detector,descriptor and matcher, create the objects with the parameters taken from respective ini files
    int assignMethod(std::string str);//given the name of the method it returns the id specified by defines on the top.


};

#endif // FEATURESELECTOR_H
