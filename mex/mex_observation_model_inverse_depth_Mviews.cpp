#include "mex.h"
#include "../src/PwgOptimiser.h"

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "triangulate requires four input arguments.");
    } else if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "triangulate requires one output argument.");
    }
    
    /* create a pointer to the real data in the input matrix  */
    const double *x = mxGetPr(prhs[0]); /* (N+6)x1 */
    double *p = mxGetPr(prhs[1]); /* (2xN)   */
    size_t ncol = mxGetN(prhs[1]);
    double *ptr = mxGetPr(prhs[2]); /* inverse depth index */
    int i = ptr[0];
    ptr = mxGetPr(prhs[3]); /* camera index */
    int c = ptr[0];
    plhs[0] = mxCreateDoubleMatrix(2*ncol, 1, mxREAL);
    double *z = mxGetPr(plhs[0]);
    
    /****************************/

    //PwgOptimiser object;
    //object.observation_model_inverse_depth(z, x, p, N);
    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object ; // pointer initialisation
    Object = new PwgOptimiser ( ) ; // pointer initialisation
    Object->observation_model_inverse_depth_Mviews(z, x, p, i-1, c-1);

}