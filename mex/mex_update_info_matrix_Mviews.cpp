#include "mex.h"
#include <iostream> // std::cout
#include <string>

#include "../src/PwgOptimiser.cpp"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
    /* create a pointer to the real data in the input matrix  */
    int nfields = mxGetNumberOfFields(prhs[0]); // number of fields in each constraint
    mwSize NStructElems = mxGetNumberOfElements(prhs[0]); // number of constraints
    const char **fnames = (const char **)mxCalloc(nfields, sizeof(*fnames)); /* pointers to field names */
    double *ptr = mxGetPr(prhs[1]); /* linearisation point */
    int ncols = ptr[0];
    ptr = mxGetPr(prhs[2]);
    int ncams = ptr[0];
    
    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object; // pointer initialisation
    Object = new PwgOptimiser; // pointer initialisation
    
    /* get constraints from Matlab to C++ and initialise a constraint */
    double *pr;
    mxArray *tmp;
    int cam, kpt;
    std::vector<double> p1(2), z(2), R(4,0.0);//, Y(49,0.0), y(7,0.0);
    std::string str1 ("cam");
    std::string str2 ("kpt");
    std::string str3 ("p1");
    std::string str4 ("z");
    std::string str5 ("R");
    std::string str6 ("y");
    std::string str7 ("Y");
    Eigen::MatrixXd Yz(7,7);
    Eigen::VectorXd yz(7);
    for (mwIndex jstruct = 0; jstruct < NStructElems; jstruct++) { /* loop the constraints */
        for (int ifield = 0; ifield < nfields; ifield++) {  /* loop the fields */
            fnames[ifield] = mxGetFieldNameByNumber(prhs[0],ifield); // get field name
            tmp = mxGetFieldByNumber(prhs[0], jstruct, ifield); // get field
            pr = (double *)mxGetData(tmp); // get the field data pointer
            if (str1.compare(fnames[ifield]) == 0) // cam
                cam = pr[0];
            if (str2.compare(fnames[ifield]) == 0) // kpt
                kpt = pr[0];
            if (str3.compare(fnames[ifield]) == 0){ // p1
                p1[0] = pr[0]; p1[1] = pr[1];}
            if (str4.compare(fnames[ifield]) == 0){ // z
                z[0] = pr[0]; z[1] = pr[1];}
            if (str5.compare(fnames[ifield]) == 0){ // R
                R[0] = pr[0]; R[3] = pr[3];}
        } // end of nfields loop
        Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz); // using pointer to object
    } // end of NStructElems loop
    
    // initialise constraints information
    Eigen::SparseMatrix<double> Y(ncols,ncols);
    Eigen::VectorXd y(ncols);
    //Object->update_info_matrix_Mviews(Y, y, xs, ncams);
    
    /* get output information vector */
    plhs[0] = mxCreateDoubleMatrix(ncols, 1, mxREAL);
    double *vec_out = mxGetPr(plhs[0]);
    for (int i=0; i<ncols; i++) {
        vec_out[i] = y.coeffRef(i);}
    
    /* get output information matrix */
    long unsigned int nzmax = Y.nonZeros();
    plhs[1] = mxCreateSparse(ncols, ncols, nzmax, mxREAL);
    pr = mxGetPr(plhs[1]);
    long unsigned int *ir = mxGetIr(plhs[1]);
    long unsigned int *jc = mxGetJc(plhs[1]);
    int index=0;
    for (int k=0; k < Y.outerSize(); ++k)
    {
        jc[k] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Y,k); it; ++it)
        {
            ir[index] = it.row();
            pr[index] = it.value();
            index++;
        }
    }
    jc[Y.outerSize()] = index;
    
    /* Free memory */
    mxFree((void *)fnames);
    delete Object; // delete class pointer
}