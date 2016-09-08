#ifndef FEATURESELECTOR_H
#define FEATURESELECTOR_H

# define    vgSLAM_KAZE    0 //det
# define    vgSLAM_FAST    1 //det
# define    vgSLAM_SIFT    2 //det & desc
# define    vgSLAM_GFTT    3 //det
# define    vgSLAM_SURF    4 //det & desc
# define    vgSLAM_BRIEF   5 //desc
# define    vgSLAM_ORB     6 //det & desc
# define    vgSLAM_BRISK   7 //det & desc
# define    vgSLAM_FREAK   8 //desc
# define    vgSLAM_FLANN   9 //matcher
# define    vgSLAM_BRUTEFORCEL2   10 //matcher
# define    vgSLAM_BRUTEFORCEL1 11 //matcher
# define    vgSLAM_BRUTEFORCEHAMMING   12//matcher
# define    vgSLAM_AKAZE    13 //det

#include "yarp/os/all.h"
#include "opencvs.h"
#include "cpps.h"


class FeatureSelector
{
    yarp::os::ResourceFinder rf;
    int argc;
    char** argv;

public:
    FeatureSelector();
    FeatureSelector(int ac, char **av);
    bool process(cv::Ptr<cv::Feature2D> &detector, cv::Ptr<cv::Feature2D> &descriptor, cv::Ptr<cv::DescriptorMatcher> &matcher);
protected:
    int  parseMap(std::string value, std::map<std::string, int> &m);
    bool checkDetector(std::string str);
    bool checkDescriptor(std::string str);
    bool checkMatcher(std::string str);
    void switcher(int detID, int descID, int matchID,
                  cv::Ptr<cv::Feature2D> &detector, cv::Ptr<cv::Feature2D> &descriptor,cv::Ptr<cv::DescriptorMatcher> &matcher);
    int assignMethod(std::string str);


};

#endif // FEATURESELECTOR_H
