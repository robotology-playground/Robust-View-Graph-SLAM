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
    Eigen::SparseMatrix<double> H(2*ncol,ncol+6);
    
    /* create and object realisation */
    //PwgOptimiser object;
    //object.observation_model_jacobian_inverse_depth(H, x, p, N);
    //std::cout << Eigen::MatrixXd(H) << "\n";
    /* initialise a PwgOptimiser object */
    //PwgOptimiser *Object; // pointer initialisation
    //Object = new PwgOptimiser ( ) ;
    //Object->observation_model_jacobian_inverse_depth(H, x, p, N);
    PwgOptimiser *Object ; // pointer initialisation
    Object = new PwgOptimiser ( ) ; // pointer initialisation
    Object->observation_model_jacobian_inverse_depth_Mviews(H, x, p, i-1, c-1);
    
    /****************************/
    
    /* get output pointers */
    long unsigned int nzmax = H.nonZeros();
    plhs[0] = mxCreateSparse(2*ncol, ncol+6, nzmax, mxREAL);
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
    
    delete Object;
}