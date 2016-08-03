#include "mex.h"
#include <iostream> // std::cout
#include <string>

#include "../src/PwgOptimiser.h"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    
    /* information vector */
    double *yin = mxGetPr(prhs[0]);
    /* sparse information matrix */
    long unsigned int ncol = mxGetN(prhs[1]);
    long unsigned int *ir = mxGetIr(prhs[1]); /* Row indexing */
    long unsigned int *jc = mxGetJc(prhs[1]); /* Column count */
    double *s = mxGetPr(prhs[1]); /* Non-zero elements */
    /* constraints structure */
    int nfields = mxGetNumberOfFields(prhs[2]); // number of fields in each constraint
    mwSize NStructElems = mxGetNumberOfElements(prhs[2]); // number of constraints
    const char **fnames = (const char **)mxCalloc(nfields, sizeof(*fnames)); /* pointers to field names */
    /* switch vector */
    double *sw = mxGetPr(prhs[3]);
    /* linearisation point */
    double *xs = mxGetPr(prhs[4]);
    long unsigned int nrow = mxGetM(prhs[4]);
    /* number of cameras */
    double *ptr = mxGetPr(prhs[5]);
    int ncams = ptr[0];
    
    /* get information matrix and vector into Eigen sparse format */
    Eigen::SparseMatrix<double> Y(ncol,ncol);
    Eigen::VectorXd y(ncol);
    for (int i=0; i<ncol; i++)
        y.coeffRef(i)=yin[i];
    std::vector<T> tripletList;
    for (int j=0; j<ncol; j++) /* Loop through columns */
        for (int i=jc[j]; i<jc[j+1]; i++) /* Loop through non-zeros in ith column */
            //std::cout << "row: " << ir[i] << ", col: " << j << ", data: " << s[i] << std::endl;
            tripletList.push_back(T(ir[i], j, s[i]));
    Y.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << Eigen::MatrixXd(Y) << std::endl;
    
    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object; // pointer initialisation
    Object = new PwgOptimiser ((int)ncams, (int)nrow-6*ncams) ;
    
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
            fnames[ifield] = mxGetFieldNameByNumber(prhs[2],ifield); // get field name
            tmp = mxGetFieldByNumber(prhs[2], jstruct, ifield); // get field
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
            if (str6.compare(fnames[ifield]) == 0) // y
                for (int ii=0; ii<7; ii++)
                    yz.coeffRef(ii)= pr[ii];
            if (str7.compare(fnames[ifield]) == 0) // Y
                for (int jj=0; jj<7; jj++)
                    for (int ii=0; ii<7; ii++)
                        Yz.coeffRef(ii,jj)= pr[jj*7+ii];
        } // end of nfields loop
        Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz, (int)sw[jstruct]); // using pointer to object
    } // end of NStructElems loop

    /* add information by constraints */
    Object->set_information_matrix_and_vector( y, Y );
    Object->constraints_addition_inverse_depth_Mviews( xs );
    Object->get_information_matrix_and_vector( y, Y );
    int sw_vec[NStructElems];
    Object->get_switch_vector( sw_vec );
    
    /* get output information vector */
    plhs[0] = mxCreateDoubleMatrix(ncol, 1, mxREAL);
    double *y_out = mxGetPr(plhs[0]);
    for (int i=0; i<ncol; i++)
        y_out[i] = y.coeffRef(i);
    
    /* get output sparse information matrix */
    long unsigned int nzmax = Y.nonZeros();
    plhs[1] = mxCreateSparse(ncol, ncol, nzmax, mxREAL);
    double *Y_out = mxGetPr(plhs[1]);
    long unsigned int *ir2 = mxGetIr(plhs[1]);
    long unsigned int *jc2 = mxGetJc(plhs[1]);
    int index=0;
    for (int k=0; k < Y.outerSize(); ++k)
    {
        jc2[k] = index;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Y,k); it; ++it)
        {
            ir2[index] = it.row();
            Y_out[index] = it.value();
            index++;
        }
    }
    jc2[Y.outerSize()] = index;
    
    /* get output switch vector */
    plhs[2] = mxCreateDoubleMatrix(1, NStructElems, mxREAL);
    double *sw_out = mxGetPr(plhs[2]);
    for (int i=0; i<NStructElems; i++)
        sw_out[i] = sw_vec[i];
    
    /* Free memory */
    mxFree((void *)fnames);
    delete Object, sw;
}