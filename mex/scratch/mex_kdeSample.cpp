#include "mex.h"
#include "../src/libs.h"
#include <algorithm>    // std::sort
#include "../src/kde.cpp"
using namespace std;

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* inputs */
    double *dim1;  /* Nx1 */
    dim1 = mxGetPr(prhs[0]);
    double *dim2;  /* Nx1 */
    dim2 = mxGetPr(prhs[1]);
    double *dim3;  /* 1   */
    dim3 = mxGetPr(prhs[2]);
    /* input dimensions */
    size_t ncols,nrows;
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if ( mxGetM(prhs[0])<mxGetN(prhs[0]) || mxGetM(prhs[1])<mxGetN(prhs[1]) ){
        mexErrMsgIdAndTxt("MATLAB:kdSample:nargin",
                "First and second inputs should be both colomn vectors.");}
    if ( mxGetN(prhs[0])>1 || mxGetN(prhs[1])>1 ){
        mexErrMsgIdAndTxt("MATLAB:kdSample:nargin",
                "First and second inputs should be both colomn vectors.");}
    /*************************************************************/
    int i, j, nd = dim3[0];
    
    KDE* kde = new KDE();
    for (i=0;i<nrows;i++){
        kde->add_data(dim1[i],dim2[i]);}
    kde->set_kernel_type(1);
    kde->set_bandwidth_opt_type(3);
    double bw1 = kde->get_bandwidth(0);
    double bw2 = kde->get_bandwidth(1);
    //cout << "# bandwidth var 1: " << bw1 << endl;
    //cout << "# bandwidth var 2: " << bw2 << endl;
    
    vector<double> p;
    srand (time(NULL));  /* initialize random seed: */
    for (i=0;i<nd;i++){
        double k = rand();
        p.push_back(k/RAND_MAX);}
     
    vector<double> samples1,samples2;
    for (i=0;i<p.size();i++){
        /////////
        //samples1.push_back(bw1*kde->NormalCDFInverse(p[i]));
        //samples2.push_back(bw2*kde->NormalCDFInverse(p[i]));
        ////////
        samples1.push_back(bw1*kde->stdnormal_inv(p[i]));
        samples2.push_back(bw2*kde->stdnormal_inv(p[i]));
    }
    /*************************************************************/
    plhs[0] = mxCreateDoubleMatrix(samples1.size(), 2, mxREAL);
    double *out0 = mxGetPr(plhs[0]);
    for (int i=0;i<samples1.size();i++){
        out0[2*i+0]=samples1[i];
        out0[2*i+1]=samples2[i];}
    /*************************************************************/
    
    delete kde;
    
}