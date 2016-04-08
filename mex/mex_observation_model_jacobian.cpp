using namespace std;

#include "../src/pwg.cpp"

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d)!=0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d)!=0.0)
#endif

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
    double *xc; /* 6x1 relative camera pose [trans rotation] */
    double *xf; /* 3xncols 3d points */
    xc = mxGetPr(prhs[0]);
    xf = mxGetPr(prhs[1]);
    
    /* get dimensions of the input matrix */
    long unsigned int nrows, ncols; /* matrix dimensions */
    ncols = mxGetN(prhs[1]);
    
    /****************************/
    
    /* put points in pwg vector format */
    int ii,jj;
    vector<pwg::point_3d> map;
    for (ii=0;ii<ncols;ii++){
        map.push_back(pwg::point_3d(xf[3*ii],xf[3*ii+1],xf[3*ii+2]));
    }
    vector<double> pose;
    for (ii=0;ii<6;ii++){
        pose.push_back(xc[ii]);
    }
    
    /****************************/
    
    SpMat H;
    pwg geometry;
    geometry.formulate_system(pose, map);
    H = geometry.observation_jacobian(ncols);
    
    /****************************/
    
    /* get output pointers */
    long unsigned int m,n;
    m = 2*ncols;
    n = 3*ncols + 6;
    //plhs[0] = mxCreateDoubleMatrix(2*ncols, 3*ncols+6, mxREAL);
    long unsigned int nzmax = H.nonZeros();
    plhs[0] = mxCreateSparse(m, n, nzmax, mxREAL);
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
}
