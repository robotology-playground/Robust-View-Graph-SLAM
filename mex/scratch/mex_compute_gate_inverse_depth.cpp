#include "mex.h"
#include <iostream> // std::cout

#include "../src/PwgOptimiser.cpp"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    /* create a pointer to the real data in the input matrix  */
    double *x = mxGetPr(prhs[0]); /* state vector */
    double *xs = mxGetPr(prhs[3]); /* linearisation point */
    int nfields = mxGetNumberOfFields(prhs[2]); // number of fields in each constraint
    mwSize NStructElems = mxGetNumberOfElements(prhs[2]); // number of constraints
    mxClassID  *classIDflags = (mxClassID *)mxCalloc(nfields, sizeof(mxClassID)); /* allocate memory  for storing classIDflags */
    const char **fnames = (const char **)mxCalloc(nfields, sizeof(*fnames)); /* pointers to field names */
    /* get the covariance matrix - used to convert from MATLAB sparse to C++ memory */
    long unsigned int ncol = mxGetN(prhs[1]);
    long unsigned int *ir = mxGetIr(prhs[1]); /* Row indexing */
    long unsigned int *jc = mxGetJc(prhs[1]); /* Column count */
    double *s  = mxGetPr(prhs[1]); /* Non-zero elements */
    
    //std::cout << "1" << std::endl;

    /* get constraints */
    /* Get input MATLAB structure fields data and convert them into C++ fields */
    double *pr;
    mxArray *tmp;
    int idx;
    std::vector<double> p1(2), p2(2), edge(2), R(4);
    
    //PwgOptimiser Object; // object initialisation
    PwgOptimiser *Object; // pointer initialisation
    Object = new PwgOptimiser; // pointer initialisation
    
    //std::cout << "2" << std::endl;
    
    for (mwIndex jstruct = 0; jstruct < NStructElems; jstruct++) { /* loop the constraints */
        for (int ifield = 0; ifield < nfields; ifield++) {  /* loop the fields */
            fnames[ifield] = mxGetFieldNameByNumber(prhs[2],ifield); // get field name
            tmp = mxGetFieldByNumber(prhs[2], jstruct, ifield); // get field
            pr = (double *)mxGetData(tmp); // get the field data pointer
            if (ifield==0){// idx
                idx = pr[0];}
            if (ifield==1){// p1
                p1[0] = pr[0]; p1[1] = pr[1];}
            if (ifield==2){// z
                p2[0] = pr[0]; p2[1] = pr[1];}
            if (ifield==3){// R
                R[0] = pr[0]; R[1] = pr[1]; R[2] = pr[2]; R[3] = pr[3];}
        } // end of nfields loop
        //Object.initialise_a_constraint(idx, p1, p2, R); // using object
        Object->initialise_a_constraint(idx, p1, p2, R); // using pointer to object
    } // end of NStructElems loop
    
    //mwSize ndim = mxGetNumberOfDimensions(prhs[2]); // number of dimension (2 for 2D array)
    //const mwSize *dims = mxGetDimensions(prhs[2]); // dim[0] and dim[1] constraint dimensions (1 x number_of_constraints)
    //long unsigned int N = dims[ndim-1];
    //std::cout << NStructElems << std::endl;
    plhs[0] = mxCreateDoubleMatrix(1, NStructElems, mxREAL);
    double *gate = mxGetPr(plhs[0]); /* output gate */
    //Object.compute_gate_inverse_depth(gate, x, ir, jc, s, xs); // using object
    Object->compute_gate_inverse_depth(gate, x, ir, jc, s, xs, ncol); // using pointer to object
    
    //std::cout << "end" << std::endl;
    
    /* Free memory */
    mxFree((void *)fnames);
    mxFree(classIDflags);
    delete Object; // delete class pointer
}