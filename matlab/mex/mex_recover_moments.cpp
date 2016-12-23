#include "mex.h"
#include <iostream> // std::cout
#include <string>

#include "RecoverMoments.h"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
    /* sparse information matrix */
    long unsigned int ncol = mxGetN(prhs[0]);
    long unsigned int *ir = mxGetIr(prhs[0]); /* Row indexing */
    long unsigned int *jc = mxGetJc(prhs[0]); /* Column count */
    double *s = mxGetPr(prhs[0]); /* Non-zero elements */
    /* information vector */
    double *yin = mxGetPr(prhs[1]);
    
    /* get information matrix and vector into Eigen format */
    Eigen::SparseMatrix<double> A(ncol,ncol);
    Eigen::VectorXd b(ncol);
    for (int i=0; i<ncol; i++)
        b.coeffRef(i)=yin[i];
    std::vector<T> tripletList;
    for (int j=0; j<ncol; j++) /* Loop through columns */
        for (int i=jc[j]; i<jc[j+1]; i++) /* Loop through non-zeros in ith column */
            //std::cout << "row: " << ir[i] << ", col: " << j << ", data: " << s[i] << std::endl;
            tripletList.push_back(T(ir[i], j, s[i]));
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    
    /****************************/
    /* initialise an object */
    RecoverMoments *Object; // pointer initialisation
    Object = new RecoverMoments(A, b); // pointer initialisation
    
    /* get output state vector */
    ncol = (Object->x).size();
    plhs[0] = mxCreateDoubleMatrix(ncol, 1, mxREAL);
    double *vec_out = mxGetPr(plhs[0]);
    for (int i=0; i<ncol; i++)
        vec_out[i] = (Object->x).coeffRef(i);
    
    /* get output sparse covariance matrix */
    ncol = (Object->P).outerSize();
    long unsigned int nzmax = (Object->P).nonZeros();
    plhs[1] = mxCreateSparse(ncol, ncol, nzmax, mxREAL);
    double *mat_out = mxGetPr(plhs[1]);
    long unsigned int *ir2 = mxGetIr(plhs[1]);
    long unsigned int *jc2 = mxGetJc(plhs[1]);
    int index=0;
    for (int k=0; k < (Object->P).outerSize(); ++k)
    {
        jc2[k] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Object->P,k); it; ++it)
        {
            ir2[index] = it.row();
            mat_out[index] = it.value();
            index++;
        }
    }
    jc2[(Object->P).outerSize()] = index;
    
}
