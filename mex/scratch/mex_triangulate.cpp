/*==========================================================
 * mex_triangulate.cpp -
 *
 *========================================================*/
/* $Revision: 1.0.0.0 $ */

#include <iostream>
#include <math.h>
#include "mex.h"
#include <limits>

using namespace std;

#include "../src/pwg.cpp"

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    /* pointers */
    double *R;  /* 3x3 rotation matrix */
    double *t;  /* 3x1 translation vector */
    double *p1; /* 2xN 2d points in first image */
    double *p2; /* 2xN 2d points in first image */
    double *xf; /* 3xN output 3d points */
    size_t nrows,ncols; /* matrix dimensions */
    
    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "triangulate requires four input arguments.");
    } else if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "triangulate requires one output argument.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    R = mxGetPr(prhs[0]);
    t = mxGetPr(prhs[1]);
    p1 = mxGetPr(prhs[2]);
    p2 = mxGetPr(prhs[3]);
    
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    
    /* get output pointers */
    plhs[0] = mxCreateDoubleMatrix(3, (mwSize)ncols, mxREAL);
    xf = mxGetPr(plhs[0]);
    
    /* an instance of pwg */
    pwg geometry;
    
    /* triangulate points */
    vector<pwg::point_3d> out = geometry.triangulate(R, t, p1, p2, (mwSize)ncols);
    
    /* assign output */
    for (int i=0; i<out.size(); i++) {
        xf[3*(i+1)-3] = out[i].x;
        xf[3*(i+1)-2] = out[i].y;
        xf[3*(i+1)-1] = out[i].z;
    }
    
}
