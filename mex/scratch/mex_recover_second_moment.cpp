
#include "../src/libs.h"
using namespace Eigen;

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* Check for proper number of arguments */
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "least-squares requires three input arguments.");}
    else if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "least-squares requires four output argument.");}
    
    /* create a pointer to the real data in the input matrix  */
    double *mat_in; /* Y */
    mat_in=mxGetPr(prhs[0]);
    
    /* get dimensions of the input matrix */
    size_t ncols;
    ncols=mxGetN(prhs[0]);
    
    /****************************/
    
    int i,j,k,idx;
    
    /****************************/
    
    // get Y as sparse //
    SparseMatrix<double, ColMajor> Y(ncols,ncols);
    std::vector<T> tripletList;
    for (i=0;i<ncols;i++){
        for (j=0;j<ncols;j++){
            idx=ncols*i+j;
            if ( abs(mat_in[idx])>0 & std::isfinite(mat_in[idx]) ){
                tripletList.push_back(T(j,i,mat_in[idx]));}}}
    Y.setFromTriplets(tripletList.begin(),tripletList.end());
    
    /* solve as  sparse */
    // Solve instead of inverse //
    SparseMatrix<double, ColMajor> I(ncols,ncols); I.setIdentity();
    SparseMatrix<double, ColMajor> sP(ncols,ncols);
    Y.makeCompressed();
    //Eigen::ConjugateGradient<SparseMatrix<double, ColMajor> > solver(Y);
    SparseQR<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver(Y);
    //SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver(Y);
    //solver.analyzePattern(Y);
    //solver.factorize(Y);
    sP=solver.solve(I);
    // convert P back to dense and output
    MatXd P=MatrixXd(sP);
    
    /****************************/
    
    plhs[0] = mxCreateDoubleMatrix(ncols, ncols, mxREAL);
    double *mat_out = mxGetPr(plhs[1]);
    for (i=0;i<ncols;i++){
        for (j=0;j<ncols;j++){
            idx=ncols*i+j;
            mat_out[idx]=P.coeffRef(j,i);}}}