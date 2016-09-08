#include "featureselector.h"

using namespace std;
using namespace cv;

map <string, int> kazeMap{{"DIFF_PM_G1", 0}, {"DIFF_PM_G2", 1}, {"DIFF_WEICKERT", 2}, {"DIFF_CHARBONNIER", 3}},
                  akazeMap{{"DESCRIPTOR_KAZE", 3}, {"DESCRIPTOR_KAZE_UPRIGHT", 2}, {"DESCRIPTOR_MLDB", 4}, {"DESCRIPTOR_MLDB_UPRIGHT", 5}},
                  fastMap{{"TYPE_5_8", 0}, {"TYPE_7_12", 1}, {"TYPE_9_16", 2}},
                  orbMap{{"HARRIS_SCORE", 0}, {"FAST_SCORE", 1}};

int FeatureSelector::parseMap(string value, map<string,int> &m){
    map <string, int>::const_iterator iValue = m.find(value);
            if (iValue  == m.end())
                throw runtime_error("Check the .ini file, wrong enum");
            return iValue->second;

}

FeatureSelector::FeatureSelector()
{
}

FeatureSelector::FeatureSelector(int ac, char** av){
    rf.setDefaultConfigFile("../../conf/vgSLAM.ini");
    argc=ac;
    argv=av;
    rf.configure(argc,argv);
}

bool FeatureSelector::process(Ptr<Feature2D> &detector, Ptr<Feature2D> &descriptor, Ptr<DescriptorMatcher> &matcher){
    int detID=vgSLAM_KAZE, descID=vgSLAM_SIFT, matchID=vgSLAM_FLANN; // default

    if(checkDetector(rf.find("Detector").asString()))
        detID=assignMethod(rf.find("Detector").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default detector:KAZE"<<endl;
    if(checkDescriptor(rf.find("Descriptor").asString()))
        descID=assignMethod(rf.find("Descriptor").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default descriptor:SIFT"<<endl;
    if(checkMatcher(rf.find("Matcher").asString()))
        matchID=assignMethod(rf.find("Matcher").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default matcher:FLANN"<<endl;

    switcher(detID,descID,matchID,detector,descriptor,matcher);


    return false;
}

bool FeatureSelector::checkDetector(string str){
    if(str.compare("KAZE")!=0 && str.compare("AKAZE")!=0  && str.compare("FAST")!=0 && str.compare("SIFT")!=0
            && str.compare("GFTT")!=0 && str.compare("SURF")!=0
            && str.compare("ORB")!=0 && str.compare("BRISK")!=0){
        cout<<"Wrong detector argument, available options are: KAZE, FAST, SIFT, GFTT, SURF, ORB, BRISK, AKAZE"<<endl;
        return false;
    }
    else
        return true;
}
bool FeatureSelector::checkDescriptor(string str){
    if(str.compare("SIFT")!=0 && str.compare("SURF")!=0 && str.compare("BRIEF")!=0
            && str.compare("ORB")!=0 && str.compare("BRISK")!=0 && str.compare("FREAK")!=0){
        cout<<"Wrong descriptor argument, available options are: SIFT, SURF, BRIEF, ORB, BRISK, FREAK"<<endl;
        return false;
    }
    else
        return true;
}
bool FeatureSelector::checkMatcher(string str){
    if(str.compare("FLANN")!=0 && str.compare("BRUTEFORCEL2")!=0
            && str.compare("BRUTEFORCEL1")!=0 && str.compare("BRUTEFORCEHAMMING")!=0){
        cout<<"Wrong matcher argument, available options are: FLANN, BRUTEFORCEL2, BRUTEFORCEL1, BRUTEFORCEHAMMING"<<endl;
        return false;
    }
    else
        return true;
}

void FeatureSelector::switcher(int detID, int descID, int matchID, Ptr<Feature2D> &detector,
                               Ptr<Feature2D> &descriptor,Ptr<DescriptorMatcher> &matcher){

    yarp::os::ResourceFinder rfdet, rfdesc, rfmatch;
    // detector switch
    switch (detID) {
    case vgSLAM_KAZE : {
        rfdet.setDefaultConfigFile("../../conf/KAZE.ini");
        rfdet.configure(argc, argv);
        detector=KAZE::create(rfdet.find("extended").asBool(),rfdet.find("upright").asBool(),
                                    (float) rfdet.find("threshold").asDouble(), rfdet.find("nOctaves").asInt(),
                              rfdet.find("nOctaveLayers").asInt(), parseMap((string)rfdet.find("diffusivity").asString(),kazeMap));
        break;
    }
    case vgSLAM_FAST : {
        rfdet.setDefaultConfigFile("../../conf/FAST.ini");
        rfdet.configure(argc, argv);
        detector=FastFeatureDetector::create(rfdet.find("threshold").asInt(),rfdet.find("nonmaxSuppression").asBool(),
                                             parseMap((string)rfdet.find("type").asString(),fastMap));
        break;
    }
    case vgSLAM_SIFT : {
        rfdet.setDefaultConfigFile("../../conf/SIFT.ini");
        rfdet.configure(argc, argv);
        detector=xfeatures2d::SIFT::create(rfdet.find("nfeatures").asInt(),rfdet.find("nOctaveLayers").asInt(),
                                                              rfdet.find("contrastThreshold").asDouble(),
                                                              rfdet.find("edgeThreshold").asDouble(),rfdet.find("sigma").asDouble());
        break;
    }
    case vgSLAM_GFTT: {
        rfdet.setDefaultConfigFile("../../conf/GFTT.ini");
        rfdet.configure(argc, argv);
        detector=GFTTDetector::create(rfdet.find("maxCorners").asInt(),rfdet.find("qualityLevel").asDouble(),
                                                    rfdet.find("minDistance").asDouble(),rfdet.find("blockSize").asInt(),
                                                    rfdet.find("useHarrisDetector").asBool(),rfdet.find("k").asDouble());


        break;
    }
    case vgSLAM_SURF : {
        rfdet.setDefaultConfigFile("../../conf/SURF.ini");
        rfdet.configure(argc, argv);
        detector=xfeatures2d::SURF::create(rfdet.find("hessianThreshold").asDouble(),rfdet.find("nOctaves").asInt(),
                                                              rfdet.find("inOctaveLayers").asInt(),rfdet.find("extended").asBool(),
                                                              rfdet.find("upright").asBool());
        break;
    }
    case vgSLAM_ORB : {
        rfdet.setDefaultConfigFile("../../conf/ORB.ini");
        rfdet.configure(argc, argv);
        detector=ORB::create(rfdet.find("nfeatures").asInt(),rfdet.find("scaleFactor").asDouble(),
                                 rfdet.find("nlevels").asInt(),rfdet.find("edgeThreshold").asInt(),
                                 rfdet.find("firstLevel").asInt(),rfdet.find("WTA_K").asInt(),
                                 parseMap((string)rfdet.find("scoreType").asString(),orbMap),rfdet.find("patchSize").asInt(),
                                 rfdet.find("fastThreshold").asInt());
        break;
    }
    case vgSLAM_BRISK : {
        rfdet.setDefaultConfigFile("../../conf/BRISK.ini");
        rfdet.configure(argc, argv);
        detector=BRISK::create(rfdet.find("thresh").asInt(),rfdet.find("octaves").asInt(),rfdet.find("patternScale").asDouble());
        break;
    }
    case vgSLAM_AKAZE : {
        rfdet.setDefaultConfigFile("../../conf/AKAZE.ini");

        rfdet.configure(argc, argv);
        detector=AKAZE::create(parseMap((string)rfdet.find("descriptor_type").asString(),akazeMap),rfdet.find("descriptor_size").asInt(),
                               rfdet.find("descriptor_channels").asInt(),
                               rfdet.find("threshold").asDouble(),rfdet.find("nOctaves").asInt(),rfdet.find("nOctaveLayers").asInt(),
                               parseMap((string)rfdet.find("diffusivity").asString(),kazeMap));
        break;
    }
    default:
        break;
    }
    // descriptor switch
    switch (descID) {
    case vgSLAM_SIFT : {
        rfdesc.setDefaultConfigFile("../../conf/SIFT.ini");
        rfdesc.configure(argc, argv);
        descriptor=xfeatures2d::SIFT::create(rfdesc.find("nfeatures").asInt(),rfdesc.find("nOctaveLayers").asInt(),
                                                              rfdesc.find("contrastThreshold").asDouble(),
                                                              rfdesc.find("edgeThreshold").asDouble(),rfdesc.find("sigma").asDouble());
        break;
    }
    case vgSLAM_SURF : {
        rfdesc.setDefaultConfigFile("../../conf/SURF.ini");
        rfdesc.configure(argc, argv);
        descriptor=xfeatures2d::SURF::create(rfdesc.find("hessianThreshold").asDouble(),rfdesc.find("nOctaves").asInt(),
                                                              rfdesc.find("inOctaveLayers").asInt(),rfdesc.find("extended").asBool(),
                                                              rfdesc.find("upright").asBool());
        break;
    }
    case vgSLAM_BRIEF : {
        rfdesc.setDefaultConfigFile("../../conf/BRIEF.ini");
        rfdesc.configure(argc, argv);
        descriptor=xfeatures2d::BriefDescriptorExtractor::create(rfdesc.find("bytes").asInt(),
                                                                  rfdesc.find("use_orientation").asBool());
        break;
    }
    case vgSLAM_ORB : {
        rfdesc.setDefaultConfigFile("../../conf/ORB.ini");
        rfdesc.configure(argc, argv);
        descriptor=ORB::create(rfdesc.find("nfeatures").asInt(),rfdesc.find("scaleFactor").asDouble(),
                                 rfdesc.find("nlevels").asInt(),rfdesc.find("edgeThreshold").asInt(),
                                 rfdesc.find("firstLevel").asInt(),rfdesc.find("WTA_K").asInt(),
                                 parseMap((string)rfdet.find("scoreType").asString(),orbMap),rfdesc.find("patchSize").asInt(),
                                 rfdesc.find("fastThreshold").asInt());
        break;
    }
    case vgSLAM_BRISK : {
        rfdesc.setDefaultConfigFile("../../conf/BRISK.ini");
        rfdesc.configure(argc, argv);
        descriptor=BRISK::create(rfdesc.find("thresh").asInt(),rfdesc.find("octaves").asInt(),rfdesc.find("patternScale").asDouble());
        break;
    }
    case vgSLAM_FREAK : {
        rfdesc.setDefaultConfigFile("../../conf/FREAK.ini");
        rfdesc.configure(argc, argv);
        descriptor=xfeatures2d::FREAK::create(rfdesc.find("orientationNormalized").asBool(),rfdesc.find("scaleNormalized").asBool(),
                                              rfdesc.find("patternScale").asDouble(),rfdesc.find("nOctaves").asInt());
        break;
    }
    default:
        break;
    }
    // matcher switch
    switch (matchID) {
    case vgSLAM_FLANN:{
        matcher=DescriptorMatcher::create("FlannBased");
        break;
    }
    case vgSLAM_BRUTEFORCEL2:{
        matcher=DescriptorMatcher::create("BruteForce");
        break;
    }
    case vgSLAM_BRUTEFORCEL1:{
        matcher=DescriptorMatcher::create("BruteForce-L1");
        break;
    }
    case vgSLAM_BRUTEFORCEHAMMING:{
        matcher=DescriptorMatcher::create("BruteForce-Hamming");
        break;
    }
    default:
        break;
    }

}


int FeatureSelector::assignMethod(string str){
    if (str.compare("KAZE")==0)
        return vgSLAM_KAZE;
    if (str.compare("FAST")==0)
        return vgSLAM_FAST;
    if (str.compare("SIFT")==0)
        return vgSLAM_SIFT;
    if (str.compare("GFTT")==0)
        return vgSLAM_GFTT;
    if (str.compare("SURF")==0)
        return vgSLAM_SURF;
    if (str.compare("BRIEF")==0)
        return vgSLAM_BRIEF;
    if (str.compare("ORB")==0)
        return vgSLAM_ORB;
    if (str.compare("BRISK")==0)
        return vgSLAM_BRISK;
    if (str.compare("FREAK")==0)
        return vgSLAM_FREAK;
    if (str.compare("FLANN")==0)
        return vgSLAM_FLANN;
    if (str.compare("BRUTEFORCEL1")==0)
        return vgSLAM_BRUTEFORCEL1;
    if (str.compare("BRUTEFORCEL2")==0)
        return vgSLAM_BRUTEFORCEL2;
    if (str.compare("BRUTEFORCEHAMMING")==0)
        return vgSLAM_BRUTEFORCEHAMMING;
    if (str.compare("AKAZE")==0)
        return vgSLAM_AKAZE;
    return -1;
}
