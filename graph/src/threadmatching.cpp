/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadmatching.h"
#include <iCub/iKin/iKinFwd.h>
#include <yarp/math/Math.h>
#include <yarp/sig/Vector.h>

using namespace iCub::ctrl;
using namespace iCub::iKin;
using namespace yarp::math;
using namespace yarp::sig;

using namespace cv;

ThreadMatching::ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn1, vgSLAMBuffer<SlamType> &bufferIn2, vgSLAMBuffer<SlamType> &bufferOut,
                               cv::Ptr<cv::DescriptorMatcher> _matcher): vgSLAMThread(bufferIn1, bufferOut),matcher(_matcher){
    first = second = NULL;
    ThreadMatching::bufferIn2 = &bufferIn2;
    ThreadMatching::bufferIn1 = &bufferIn1;

}
ThreadMatching::~ThreadMatching(){
    sampler.save("./TimeStamp.log");
    if(first) {
        first->free();
        delete first;
    }

    if(second) {
        second->free();
        delete second;
    }
}

void ThreadMatching::interrupt() {
    bufferIn2->interrupt();
    vgSLAMThread::interrupt();
}


boost::dynamic_bitset<> ThreadMatching::remove_points_at_infinity(std::vector<Point2f>  matched1, std::vector<Point2f>  matched2, double thres) {
    /* matches should move at least pixel_disparity pixels
     * the largest the pixel_disparity threshold is, the more stable the results
     * and the estimates, but the less dense the result is and the shorter the
     * observed range is.
     */
    float diff1, diff2;
    boost::dynamic_bitset<> vis(matched1.size()); // all 0's by default
    for (int i=0; i<matched1.size(); i++){
        diff1 = matched2[i].x - matched1[i].x;
        diff2 = matched2[i].y - matched1[i].y;
        if ( fabs(diff1) > thres | fabs(diff2) > thres )
            vis[i] = 1;
    }
    return vis;
}

//it computes the projection matrix using the angles of the encoders read from the ports.
void ThreadMatching::getKinTransformationsToRoot(Mat &ProjectionMatrix,SlamType* data){
    std::string str="";
    if(data->right){
        str="right";}
    else{
        str="left";}
    iCubEye eye(str);
    Vector q0,qf,qhat,xf,xhat;

    Matrix trasf;

    iKinChain *chain;

    chain=eye.asChain();

    chain->setAllConstraints(false);//Ugo suggest to drop the kinematics contraint because we are not connected to a robot for now

    q0=chain->getAng();
    std::cout << "Unblocking the torso joints... "<<std::endl;
    chain->releaseLink(0);
    chain->releaseLink(1);
    chain->releaseLink(2);

    std::cout << chain->getDOF() << " DOFs available" << std::endl;

    qf.resize(chain->getDOF());
    double version,vergence;

    version=data->anglesHead->at(4)*CTRL_DEG2RAD;
    vergence=data->anglesHead->at(5)*CTRL_DEG2RAD;

    qf[0]=data->anglesTorso->at(2)*CTRL_DEG2RAD;//the torso angles are inverted
    qf[1]=data->anglesTorso->at(1)*CTRL_DEG2RAD;
    qf[2]=data->anglesTorso->at(0)*CTRL_DEG2RAD;
    qf[3]=data->anglesHead->at(0)*CTRL_DEG2RAD;
    qf[4]=data->anglesHead->at(1)*CTRL_DEG2RAD;
    qf[5]=data->anglesHead->at(2)*CTRL_DEG2RAD;
    qf[6]=data->anglesHead->at(3)*CTRL_DEG2RAD;//tilt, the eye have two DOF, tilt and pan
    qf[7]=version - vergence/2;//pan

    trasf=chain->getH(qf);

    ProjectionMatrix = (cv::Mat_<double>(4,4) << trasf(0,0), trasf(0,1), trasf(0,2), trasf(0,3),
                        trasf(1,0), trasf(1,1), trasf(1,2), trasf(1,3),
                        trasf(2,0), trasf(2,1), trasf(2,2), trasf(2,3),
                        trasf(3,0), trasf(3,1), trasf(3,2), trasf(3,3));




}

cv::Mat ThreadMatching::getProjMat(MatchesVector &matches, SlamType* data1, SlamType *data2){
    Mat Proj(4,4,CV_64F);
    double max_dist = 0;
    double min_dist = 100;
    double mil_dist = 5;
    //TESTED it works, correct order.
//    if(data1->right){
//        yError()<<"data1:RIGHT";}
//    else{
//        yError()<<"data1:LEFT";}

//    if(data2->right){
//        yError()<<"data2:RIGHT";}
//    else{
//        yError()<<"data2:LEFT";}

    //-- Quick calculation of max and min distances between keypoints
    // http://docs.opencv.org/3.1.0/d5/d6f/tutorial_feature_flann_matcher.html#gsc.tab=0
    for(int i = 0;i<matches.size();i++){
        double dist = matches[i].distance;
        if(dist<min_dist) min_dist=dist;
        if(dist>max_dist) max_dist=dist;
    }
    //cout<<"-- Max dist : " <<max_dist<<"\n";
    //cout<<"-- Min dist : " <<min_dist<<"\n";

    //-- Get only "good" matches (i.e. whose distance is less than 2*min_dist,
    //-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
    //-- small)
    //-- PS.- radiusMatch can also be used here.
    MatchesVector good_matches;
    for(int i=0;i<matches.size();i++)
        if(matches[i].distance<=max(mil_dist*min_dist,0.02))
            good_matches.push_back(matches[i]);

    //-- Collect the matching key-points
    std::vector<Point2f> kp1,kp2,matched1, matched2;

    KeyPoint::convert(*data1->feature,kp1);
    KeyPoint::convert(*data2->feature,kp2);

    for(int i=0;i<good_matches.size();i++){ // using VLMATCH
        matched1.push_back(kp1[good_matches[i].queryIdx]);
        matched2.push_back(kp2[good_matches[i].trainIdx]);
    }

    //-- Remove small tracks
    boost::dynamic_bitset<> vis;
    vis = remove_points_at_infinity(matched1, matched2, 2);
    std::vector<Point2f>  vismatched1, vismatched2;
    for(int i = 0; i < matches.size(); i++)
        if(vis[i]){
            vismatched1.push_back(matched1[i]);
            vismatched2.push_back(matched2[i]);
        }
    yInfo() << "ThreadMatching: visible matches= " << (int)vis.count();
    if (vis.count()>100){
        //-- http://nghiaho.com/?p=1675
        //-- https://avisingh599.github.io/vision/monocular-vo/
        Mat K1 = (Mat_<double>(3,3) << 235.1162, 0, 155.6680, 0, 235.7933, 122.0000, 0, 0, 1);
        Mat K2 = (Mat_<double>(3,3) << 234.2173, 0, 149.1942, 0, 234.8432, 124.5134, 0, 0, 1);
        Mat k1 = (Mat_<double>(1,5) << -0.4323, 0.1952, -0.0003,  0.0018, 0);
        Mat k2 = (Mat_<double>(1,5) << -0.4309, 0.1880,  0.0007, -0.0008, 0);

        //-- Normalize
        undistortPoints(vismatched1, vismatched1, K1, k1, noArray(), K1);
        undistortPoints(vismatched2, vismatched2, K2, k2, noArray(), K2);


        //-- Essential matrix with RANSAC
        Mat E, R, t, mask;
        E = findEssentialMat(vismatched1, vismatched2, 1.0, Point2d(0,0), RANSAC, 0.999, 0.0001, mask);
        //correctMatches(E, vismatched1, vismatched2, vismatched1, vismatched2);
        recoverPose(E, vismatched1, vismatched2, R, t, 1.0, Point2d(0,0), noArray());
        Mat a;
        Rodrigues(R,a);
        transpose(a,a);
        Mat m=(cv::Mat_<double>(1,4) << 0.0, 0.0, 0.0, 1.0);
        //std::cout<<"ThreadMatching: translation "<< t <<"angle"<< a*180/M_PI<<std::endl;
        hconcat(R,t,Proj);
        vconcat(Proj,m,Proj);
        yInfo()<<"ThreadMatching::I'm using vision for the projection matrix";
        std::cout<</*R<<std::endl<<t<<std::endl<<*/Proj<<std::endl;
        return Proj;
    }
        //in case we have only few good matches we compute Proj using kinematics.
    else{
        Mat Proj1(4,4,CV_64F),Proj2(4,4,CV_64F);
        yInfo()<<"ThreadMatching::I'm using kinematics for the projection matrix";
        //kinematics
        getKinTransformationsToRoot(Proj1,data1);
        getKinTransformationsToRoot(Proj2,data2);
        Proj=Proj1.inv()*Proj2;//to check, I think it is correct.
        std::cout<<Proj<<std::endl;
        return Proj;
    }
    return Proj;
}

void ThreadMatching::run (){
    bool change=false;
    //We continue to read, process and write data while it is not interrupted and there is data to process.
    //The interrupted flag is setted to true in the function interrupt() of the parent class "vgSLAMThread".
    //interrupt() is called by thread.close() when we close the module.
    while(!interrupted) {
        //The thickness now is 3, we have bundles of 3 frames, at first loop we compute 1-2,2-3,1-3
        //then in the others we discard the first data(1), we take a new one(3) and we compute only 1-3 and 2-3
        //because 1-2 was already computed in the previous step(in the previous was 2-3).
        yInfo()<<"ThreadMatching:Reading buffer";
        SlamType* data = new SlamType();
        MatchesVector mv;
        //Every loop it reads a new frame, one time from the L buffer, one time from R buffer alternatively
        if (change){
            change=false;
            if(!bufferIn2->read(*data)) {
                delete data;
                continue;
            }
        }
        else {
            change=true;
            if(!bufferIn1->read(*data)) {
                delete data;
                continue;
            }
        }
        sampler.add(data->stamp->getTime());
        if(!first) {
            first = data;
            continue;
        }

        if(!second) {
            second = data;
            data->matching=new MatchesVector;
            matcher->match(*first->descriptor,*second->descriptor,mv,cv::noArray());
            yDebug()<<"ThreadMatching first->second: size matching"<<mv.size();
            getProjMat(mv,first,second);
            continue;
        }

        // processing
        yDebug()<<"ThreadMatching: first:"<<first->stamp->getCount()<<"second:"<<second->stamp->getCount()<<"third:"<<data->stamp->getCount();//tested ok
        // ...
        matcher->match(*second->descriptor,*data->descriptor,mv,cv::noArray());
        getProjMat(mv,second,data);
        yDebug()<<"ThreadMatching second->third: size matching"<<mv.size();
        matcher->match(*first->descriptor,*data->descriptor,mv,cv::noArray());
        getProjMat(mv,first,data);
        yDebug()<<"ThreadMatching first->third: size matching"<<mv.size();


        //substitution
        first->free();
        delete first;
        first = second;
        second = data;
        countProcessed++;

    } //end while
}
