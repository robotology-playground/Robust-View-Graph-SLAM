#include <iostream>
#include <string.h>

#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include "mexopencv.hpp"
//#include "mex.h"
#include "../../cpp/icub/src/image.h"

using namespace std;
using namespace cv;

void mexFunction (int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    // Check the number of arguments
    //if (nrhs<3 || ((nrhs%2)!=1) || nlhs>1)
    //    mexErrMsgIdAndTxt("Wrong number of arguments");
    
    // Argument vector
    vector<MxArray> rhs(prhs,prhs+nrhs);
    Mat img1(rhs[0].toMat());
    Mat img2(rhs[1].toMat());
    
    // Option processing
    //Mat newCameraMatrix;
    //for (int i=3; i<nrhs; i+=2) {
    //    string key = rhs[i].toString();
    //    if (key=="NewCameraMatrix")
    //        newCameraMatrix = rhs[i+1].toMat(CV_32F);
    //    else
    //        mexErrMsgIdAndTxt("mexopencv:error","Unrecognized option");
    //}
    
    // Process
    image object;
    object.features(img1);  // process the first image
    object.features(img2);  // process the second image
    
}