/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadmatching.h"

using namespace cv;

ThreadMatching::ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,
                               cv::Ptr<cv::DescriptorMatcher> _matcher): vgSLAMThread(bufferIn, bufferOut),matcher(_matcher){
    first = second = NULL;

}
ThreadMatching::~ThreadMatching(){
    if(first) {
        first->free();
        delete first;
    }

    if(second) {
        second->free();
        delete second;
    }
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

cv::Mat ThreadMatching::getProjMat(MatchesVector &matches, SlamType* data1, SlamType *data2){
    Mat Proj(4,3,CV_64F);
    double max_dist = 0;
    double min_dist = 100;
    double mil_dist = 5;
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
        //std::cout<<"ThreadMatching: translation "<< t <<"angle"<< a*180/M_PI<<std::endl;
        hconcat(R,t,Proj);
        std::cout<<R<<std::endl<<t<<std::endl<<Proj<<std::endl;
    }
    return Proj;
}

void ThreadMatching::run (){

    while(!interrupted) {
        yInfo()<<"ThreadMatching:Reading buffer";
        SlamType* data = new SlamType();
        MatchesVector mv;
        if(!bufferIn->read(*data)) {
            delete data;
            continue;
        }

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
