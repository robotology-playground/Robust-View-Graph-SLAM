/*==========================================================
 * mex_least_squares.cpp -
 *
 *========================================================*/
/* $Revision: 1.0.0.0 $ */

using namespace std;

#include "pwg.cpp"

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    /* pointers */
    int i;
    double *R;  /* 3x3 rotation matrix - input */
    double *t;  /* 3x1 translation vector - input  */
    double *p1; /* 2xN 2d points in first image - input  */
    double *p2; /* 2xN 2d points in second image - input  */
    double *X;  /* 2xN 3d points - input */
    double *xf; /* 3xN output 3d points - output  */
    double *xc; /* 6x1 output 3d points - output */
    double *trc; /* 6x1 trace - output */
    double *sw; /* 2*Nx1 switch - output */
    size_t nrows,ncols; /* matrix dimensions */
    
    /* Check for proper number of arguments */
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "triangulate requires four input arguments.");
    } else if (nlhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "triangulate requires one output argument.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    R = mxGetPr(prhs[0]);
    t = mxGetPr(prhs[1]);
    p1 = mxGetPr(prhs[2]);
    p2 = mxGetPr(prhs[3]);
    X = mxGetPr(prhs[4]);
    
    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[2]);
    
    /* get output pointers */
    plhs[0] = mxCreateDoubleMatrix(6, 1, mxREAL);
    xc = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(3, (mwSize)ncols, mxREAL);
    xf = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(6, 1, mxREAL);
    trc = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1, 2*(mwSize)ncols, mxREAL);
    sw = mxGetPr(plhs[3]);
    
    /****************************/
    
    /* put points in pwg vector format */
    vector<pwg::point_3d> map;
    for (i=0;i<ncols;i++){
        map.push_back(pwg::point_3d(X[3*i-3],X[3*i-2],X[3*i-1]));
    }
    vector<double> pose;
    for (i=0;i<3;i++){
        pose.push_back(t[i]);
    }
    
    /****************************/
    
    /* an instance of pwg */
    pwg geometry;
    
    /* build the sparse system of equations */
    geometry.formulate_system(pose, map);
    
    /* look at the matrices */
    //libcommon obj;
    //obj.print_matrix_sparse(H);
    //obj.print_matrix_sparse(H);
    
    /****************************/
    
    /* assign output */
    for (int i=0; i<map.size(); i++) {
        xf[3*(i+1)-3] = map[i].x;
        xf[3*(i+1)-2] = map[i].y;
        xf[3*(i+1)-1] = map[i].z;
    }
    
}
