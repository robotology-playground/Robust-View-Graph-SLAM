using namespace std;

#include "../src/pwg.cpp"

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "least-squares requires three input arguments.");} 
    else if (nlhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "least-squares requires four output argument.");}
    
    /* create a pointer to the real data in the input matrix  */
    double *pix;/* 2xncols 2d points in second image */
    double *xc; /* 6x1 relative camera pose [trans rotation] */
    double *xf; /* 3xncols 3d points */
    pix = mxGetPr(prhs[0]);
    xc = mxGetPr(prhs[1]);
    xf = mxGetPr(prhs[2]);
    
    /* get dimensions of the input matrix */
    size_t nrows, ncols; /* matrix dimensions */
    ncols = mxGetN(prhs[2]);
    
    /****************************/
    
    /* put points in pwg vector format */
    size_t m,n;
    m = 2*ncols;
    n = 3*ncols + 6;
    int i;
    VecXd obs(m);
    for (i=0;i<ncols;i++){
        obs.coeffRef(2*i+0) = pix[2*i+0];
        obs.coeffRef(2*i+1) = pix[2*i+1];}
    vector<pwg::point_3d> map;
    for (i=0;i<ncols;i++){
        map.push_back(pwg::point_3d(xf[3*i+0],xf[3*i+1],xf[3*i+2]));}
    vector<double> pose;
    for (i=0;i<6;i++){
        pose.push_back(xc[i]);}
    
    /****************************/
    
    VecXd zhat(m);
    SpMat H;
    
    /* an instance of pwg */
    pwg geometry;
    
    /* build the sparse system of equations */
    geometry.formulate_system(pose, map);
    
    /* compute measurements and jacobian */
    zhat = geometry.observation_model(ncols);
    H = geometry.observation_jacobian(ncols);
    
    /* adding information */
    geometry.constraints_addition(obs, zhat, H);
    
    /* removing information */
    geometry.constraints_removal(obs, zhat, H);
    
    
    /****************************/
    
    /* get output pointers */
    plhs[0] = mxCreateDoubleMatrix(6, 1, mxREAL);
    double *xc_hat = mxGetPr(plhs[0]);
    for (i=0; i<6; i++) {
        xc_hat[i] = geometry.filter.x.coeffRef(3*ncols+i);}
    
    plhs[1] = mxCreateDoubleMatrix(3, ncols, mxREAL);
    double *xf_hat = mxGetPr(plhs[1]);
    for (int i=0; i<ncols; i++) {
        xf_hat[3*i+0] = geometry.filter.x.coeffRef(3*i+0);
        xf_hat[3*i+1] = geometry.filter.x.coeffRef(3*i+1);
        xf_hat[3*i+2] = geometry.filter.x.coeffRef(3*i+2);}
    
    SpMat Y = geometry.filter.Y;
    n = Y.outerSize();
    size_t nzmax = Y.nonZeros();
    plhs[2] = mxCreateSparse(n, n, nzmax, mxREAL);
    double *pr = mxGetPr(plhs[2]);
    size_t *ir = mxGetIr(plhs[2]);
    size_t *jc = mxGetJc(plhs[2]);
    int index=0;
    for (i=0; i < n; ++i){
        jc[i] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Y,i); it; ++it){
            ir[index] = it.row();
            pr[index] = it.value();
            index++;}}
    jc[n] = index;
    
    plhs[3] = mxCreateDoubleMatrix(1, 2*ncols, mxREAL);
    double *sw = mxGetPr(plhs[3]);
    for (i=0; i<2*ncols; i++) {
        sw[i] = geometry.filter.sw[i];}
    
}
