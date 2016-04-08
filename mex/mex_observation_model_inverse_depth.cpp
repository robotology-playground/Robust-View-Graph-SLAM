#include "mex.h"
#include "../src/PwgOptimiser.cpp"

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "triangulate requires four input arguments.");
    } else if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "triangulate requires one output argument.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    double *x = mxGetPr(prhs[0]); /* (N+6)x1 */
    double *p = mxGetPr(prhs[1]); /* (2xN)   */
    long unsigned int N = mxGetN(prhs[1]); /* matrix dimensions */
    plhs[0] = mxCreateDoubleMatrix(2*N, 1, mxREAL);
    double *z = mxGetPr(plhs[0]);
    
    /****************************/

    PwgOptimiser object;
    object.observation_model_inverse_depth(z, x, p, N);

}