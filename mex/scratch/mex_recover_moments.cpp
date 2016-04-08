
#include "../src/libs.h"
using namespace Eigen;

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargin",
                "least-squares requires three input arguments.");}
    else if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:triangulate:nargout",
                "least-squares requires four output argument.");}
    
    /* create a pointer to the real data in the input matrix  */
    double *vec_in; /* y */
    double *mat_in; /* Y */
    vec_in=mxGetPr(prhs[0]);
    mat_in=mxGetPr(prhs[1]);
    
    /* get dimensions of the input matrix */
    size_t ncols;
    ncols=mxGetN(prhs[1]);
    
    /****************************/
    
    int i,j,k,idx;
    
    // get y as vector //
    VectorXd y(ncols);
    for (i=0;i<ncols;i++){
        y.coeffRef(i) = vec_in[i];}
    
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
    
    /****************************/
    
//     // get Y as dense //
//     MatXd dY(ncols,ncols);
//     for (i=0;i<ncols;i++){
//         for (j=0;j<ncols;j++){
//             idx=ncols*i+j;
//             if ( abs(mat_in[idx])>0 & std::isfinite(mat_in[idx]) ){
//                 dY.coeffRef(j,i)=mat_in[idx];}}}
//     // convert Y as sparse //
//     SpMat Y=dY.sparseView();
//     std::cout << "NNZ = " << Y.nonZeros() << std::endl;
    
    /****************************/
    
    /* solve as  sparse */
    // Solve instead of inverse //
    SparseMatrix<double, ColMajor> I(ncols,ncols); I.setIdentity();
    SparseMatrix<double, ColMajor> sP(ncols,ncols);
    VectorXd x(ncols);
    
    Y.makeCompressed();
    //Eigen::ConjugateGradient<SparseMatrix<double, ColMajor> > solver(Y);
    //SparseQR<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver(Y);
    SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver(Y);
    // Compute the ordering permutation vector from the structural pattern of Y
    //solver.analyzePattern(Y);
    // Compute the numerical factorization
    //solver.factorize(Y);
    sP=solver.solve(I);
    x=solver.solve(y);
    //x=sP*y;
    
//     // Solving using Cholesky
//     // performs a Cholesky factorization of Y
//     Eigen::SimplicialCholesky<SpMat> chol(Y);
//     // use the factorization to solve for the given right hand side
//     SpMat sP = chol.solve(I);
//     VecXd x = chol.solve(y);
    
    // convert P back to dense and output
    MatXd P=MatrixXd(sP);
    
    /****************************/
    
//     /* solve as dense */
//     MatXd I = Eigen::MatrixXd::Identity(ncols,ncols);
//     MatXd P=dY.householderQr().solve(I);
//     VecXd x=P*y;
    
//     /* solve as dense per column */
//     MatXd I = Eigen::MatrixXd::Identity(ncols,ncols);
//     MatXd P(ncols,ncols);
//     for (j=0;j<ncols;j++){
//         P.col(j)=dY.householderQr().solve(I.col(j));}
//     VecXd x=P*y;
    
    /****************************/
    
    /* get output pointers */
    plhs[0] = mxCreateDoubleMatrix(ncols, 1, mxREAL);
    double *vec_out = mxGetPr(plhs[0]);
    for (i=0; i<ncols; i++) {
        vec_out[i] = x.coeffRef(i);}
    
    plhs[1] = mxCreateDoubleMatrix(ncols, ncols, mxREAL);
    double *mat_out = mxGetPr(plhs[1]);
    for (i=0;i<ncols;i++){
        for (j=0;j<ncols;j++){
            idx=ncols*i+j;
            mat_out[idx]=P.coeffRef(j,i);}}}