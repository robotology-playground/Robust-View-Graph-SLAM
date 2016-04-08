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
    
    /* create input and output data  */
    double *x = mxGetPr(prhs[0]); /* (N+6)x1 */
    double *p = mxGetPr(prhs[1]); /* (2xN)   */
    long unsigned int N = mxGetN(prhs[1]); /* number of points */
    Eigen::SparseMatrix<double> H(2*N,N+6);
    
    /* create and object realisation */
    PwgOptimiser object;
    object.observation_model_jacobian_inverse_depth(H, x, p, N);
    //std::cout << Eigen::MatrixXd(H) << "\n";
    
    /****************************/
    
    /* get output pointers */
    long unsigned int nzmax = H.nonZeros();
    plhs[0] = mxCreateSparse(2*N, N+6, nzmax, mxREAL);
    double *pr = mxGetPr(plhs[0]);
    long unsigned int *ir = mxGetIr(plhs[0]);
    long unsigned int *jc = mxGetJc(plhs[0]);
    int index=0, k;
    for (k=0; k < H.outerSize(); ++k)
    {
        jc[k] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(H,k); it; ++it)
        {
            ir[index] = it.row();
            pr[index] = it.value();
            index++;
        }
    }
    jc[H.outerSize()] = index;
    
    //delete object;
}