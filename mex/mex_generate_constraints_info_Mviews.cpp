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
    double *xs = mxGetPr(prhs[1]); /* linearisation point */
    long unsigned int nrow = mxGetM(prhs[1]);
    double *ncams = mxGetPr(prhs[2]);
    
    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object; // pointer initialisation
    Object = new PwgOptimiser ((int)ncams[0], (int)nrow-6*ncams[0]) ; // pointer initialisation
    
    /* get constraints from Matlab to C++ and initialise a constraint */
    double *pr;
    mxArray *tmp;
    int cam, kpt, sw = 0;
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
        Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz, sw); // using pointer to object
    } // end of NStructElems loop
    
    // initialise constraints information
    Object->generate_constraints_info_Mviews( xs );
    
    // pull-out constraints (Only to use private constraints in Matlab)
    std::vector<PwgOptimiser::pulled_constraint> C;
    Object->pull_constraints_Mviews(C);
    
    // output constraints in Matlab format
    /* create a 1xNStructElems struct with nfields number of fields */
    int nfields2 = 7;
    const char **fnames2 = (const char **)mxCalloc(nfields2, sizeof(*fnames2)); /* pointers to field names */
    // gather new structure field names
    for (int ifield=0; ifield<nfields; ifield++)
        fnames2[ifield] = fnames[ifield];
    fnames2[nfields+0] = "y";
    fnames2[nfields+1] = "Y";
    plhs[0] = mxCreateStructMatrix(1, NStructElems, nfields2, fnames2);
    //int ndim = mxGetNumberOfDimensions(prhs[2]);
    //const mwSize *dims = mxGetDimensions(prhs[2]);
    mxArray *fout;
    double *vec_out;
    for (mwIndex jstruct=0; jstruct<NStructElems; jstruct++) {
        for(int ifield=0; ifield<nfields2; ifield++) {
            fout = mxCreateDoubleMatrix(1, 1, mxREAL);
            if (str1.compare(fnames2[ifield]) == 0){ // cam
                fout = mxCreateDoubleMatrix(1,1,mxREAL);
                vec_out = mxGetPr(fout);
                vec_out[0] = C[jstruct].cam;
            }
            if (str2.compare(fnames2[ifield]) == 0){ // kpt
                fout = mxCreateDoubleMatrix(1,1,mxREAL);
                vec_out = mxGetPr(fout);
                vec_out[0] = C[jstruct].kpt;
            }
            if (str3.compare(fnames2[ifield]) == 0){ // p1
                fout = mxCreateDoubleMatrix(2,1,mxREAL);
                vec_out = mxGetPr(fout);
                for (int ii=0; ii<2; ii++) {
                    vec_out[ii] = C[jstruct].p1[ii];
                }
            }
            if (str4.compare(fnames2[ifield]) == 0){ // z
                fout = mxCreateDoubleMatrix(2,1,mxREAL);
                vec_out = mxGetPr(fout);
                for (int ii=0; ii<2; ii++) {
                    vec_out[ii] = C[jstruct].z[ii];
                }
            }
            if (str5.compare(fnames2[ifield]) == 0){ // R
                fout = mxCreateDoubleMatrix(2,2,mxREAL);
                vec_out = mxGetPr(fout);
                for (int ii=0; ii<4; ii++) {
                    vec_out[ii] = C[jstruct].R[ii];
                }
            }
            if (str6.compare(fnames2[ifield]) == 0){ // y
                fout = mxCreateDoubleMatrix(7,1,mxREAL);
                vec_out = mxGetPr(fout);
                for (int ii=0; ii<7; ii++) {
                    vec_out[ii] = C[jstruct].y.coeffRef(ii);
                }
            }
            if (str7.compare(fnames2[ifield]) == 0){ // Y
                fout = mxCreateDoubleMatrix(7,7,mxREAL);
                vec_out = mxGetPr(fout);
                for (int ii=0; ii<49; ii++) {
                    vec_out[ii] = C[jstruct].Y.coeffRef(ii);
                }
            }
            /* set each field in output structure */
            mxSetFieldByNumber(plhs[0], jstruct, ifield, fout);
        }
    }
    
    /* Free memory */
    mxFree((void *)fnames);
    mxFree((void *)fnames2);
    delete Object; // delete class pointer
}