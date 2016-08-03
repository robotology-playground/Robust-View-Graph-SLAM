using namespace std;

#include "../src/pwg.cpp"

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
    vector<pwg::point_3d> map;
    int i;
    for (i=0;i<ncols;i++){
        map.push_back(pwg::point_3d(xf[3*i],xf[3*i+1],xf[3*i+2]));
    }
    vector<double> pose;
    for (i=0;i<6;i++){
        pose.push_back(xc[i]);
    }
    
    /****************************/
    
    VecXd zhat;
    pwg geometry;
    geometry.formulate_system(pose, map);
    zhat = geometry.observation_model(ncols);
    
    /****************************/
    
    /* get output pointers */
    long unsigned int m,n;
    m = 2*ncols;
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *out = mxGetPr(plhs[0]);
    for (int i=0; i<zhat.size(); i++) {
        out[i] = zhat.coeffRef(i);
    }
    
}
