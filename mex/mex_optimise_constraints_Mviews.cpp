#include "mex.h"
#include <iostream> // std::cout
#include <string>

#include "../src/PwgOptimiser.h"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
    /* create a pointer to the real data in the input matrix  */
    int nfields = mxGetNumberOfFields(prhs[0]); // number of fields in each constraint
    mwSize NStructElems = mxGetNumberOfElements(prhs[0]); // number of constraints
    const char **fnames = (const char **)mxCalloc(nfields, sizeof(*fnames)); /* pointers to field names */
    double *sw = mxGetPr(prhs[1]); /* switch vector */ 
    double *xs = mxGetPr(prhs[2]); /* linearisation point */
    long unsigned int nrow = mxGetM(prhs[2]);
    double *ptr = mxGetPr(prhs[3]); /* number of cameras */
    int ncams = ptr[0];
    
    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object; // pointer initialisation
    Object = new PwgOptimiser (ncams, nrow-6*ncams) ; // pointer initialisation
    
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
    Eigen::MatrixXd yz = Eigen::MatrixXd::Zero(7,1);
    Eigen::VectorXd Yz = Eigen::MatrixXd::Zero(7,7);
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
        Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz, (int)sw[jstruct]) ; // using pointer to object
    } // end of NStructElems loop
    
    // optimise constraints information
    Object->optimise_constraints_image_inverse_depth_Mviews( xs ) ;
            
    /* get output state vector */
    int ncol = (Object->xhat).size();
    plhs[0] = mxCreateDoubleMatrix(ncol, 1, mxREAL);
    double *vec_out = mxGetPr(plhs[0]);
    for (int i=0; i<ncol; i++)
        vec_out[i] = (Object->xhat).coeffRef(i);
    
    /* get output sparse covariance matrix */
    ncol = (Object->Phat).outerSize();
    long unsigned int nzmax = (Object->Phat).nonZeros();
    plhs[1] = mxCreateSparse(ncol, ncol, nzmax, mxREAL);
    double *mat_out = mxGetPr(plhs[1]);
    long unsigned int *ir2 = mxGetIr(plhs[1]);
    long unsigned int *jc2 = mxGetJc(plhs[1]);
    int index=0;
    for (int k=0; k < (Object->Phat).outerSize(); ++k)
    {
        jc2[k] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Object->Phat,k); it; ++it)
        {
            ir2[index] = it.row();
            mat_out[index] = it.value();
            index++;
        }
    }
    jc2[(Object->Phat).outerSize()] = index;
    
    /* Free memory */
    mxFree((void *)fnames);
    delete Object; // delete class pointer
    //mxFree(tmp);
}