#include "mex.h"
#include <iostream> // std::cout

#include "../src/GraphOptimiser.cpp"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
    /* create a pointer to the real data in the input matrix  */
    double *x = mxGetPr(prhs[0]); /* state vector */
    double *xs = mxGetPr(prhs[3]); /* linearisation point */
    int nfields = mxGetNumberOfFields(prhs[2]);             // number of fields in each constraint
    mwSize NStructElems = mxGetNumberOfElements(prhs[2]);   // number of constraints
    mxClassID  *classIDflags = (mxClassID *)mxCalloc(nfields, sizeof(mxClassID)); /* allocate memory  for storing classIDflags */
    const char **fnames = (const char **)mxCalloc(nfields, sizeof(*fnames)); /* pointers to field names */
    
    /* get the covariance matrix - used to convert from MATLAB sparse to C++ memory */
    double *s  = mxGetPr(prhs[1]); /* Non-zero elements */  // incase P is sparse or dense
    long unsigned int ncols = mxGetN(prhs[1]);              // incase P is sparse or dense
    
    /* get constraints */
    /* Get input MATLAB structure fields data and convert them into C++ fields */
    double *pr;
    mxArray *tmp;
    std::vector<double> edge(2), z(6), R(36);
    
    GraphOptimiser *Object; // pointer initialisation
    Object = new GraphOptimiser; // pointer initialisation
    
    for (mwIndex jstruct = 0; jstruct < NStructElems; jstruct++) { /* loop the constraints */
        for (int ifield = 0; ifield < nfields; ifield++) {  /* loop the fields */
            fnames[ifield] = mxGetFieldNameByNumber(prhs[2],ifield); // get field name
            tmp = mxGetFieldByNumber(prhs[2], jstruct, ifield); // get field
            pr = (double *)mxGetData(tmp); // get the field data pointer
            if (ifield==0){// edge
                edge[0] = pr[0]; 
                edge[1] = pr[1];
            }
            if (ifield==1){// z
                z[0] = pr[0]; 
                z[1] = pr[1];
                z[2] = pr[2];
                z[3] = pr[3];
                z[4] = pr[4];
                z[5] = pr[5];
            }
            if (ifield==2){// R
                R[0] = pr[0]; R[7] = pr[7]; R[14] = pr[14]; 
                R[21] = pr[21]; R[28] = pr[28]; R[35] = pr[35];
            }
        } // end of nfields loop
        Object->initialise_a_constraint(edge, z, R); // using pointer to object
    } // end of NStructElems loop
    
    plhs[0] = mxCreateDoubleMatrix(1, NStructElems, mxREAL);
    double *gate = mxGetPr(plhs[0]); /* output gate */
    Object->compute_gate(gate, x, s, xs, ncols); // using pointer to object
    
    /* Free memory */
    mxFree((void *)fnames);
    mxFree(classIDflags);
    delete Object; // delete class pointer
}