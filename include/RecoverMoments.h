#ifndef RECOVERMOMENTS_H
#define RECOVERMOMENTS_H

#include "matlab/mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream> // std::cout
#include <time.h>
#include <unistd.h> // sleep command

#include "Eigen/Sparse"
typedef Eigen::Triplet<double> T;
#include "matlab/lapack.h"
#include "matlab/blas.h"
#include "cholmod.h"
//#include "UFconfig.h" // GPStuff
#include "SuiteSparse_config.h" // SuiteSparse
#ifndef NPARTITION
#include <cholmod_partition.h>
#endif
#ifndef NSUPERNODAL
#include <cholmod_supernodal.h>
#endif

//#define Int UF_long // GPStuff
#define Int SuiteSparse_long // SuiteSparse
#define TRUE 1
#define FALSE 0
#define CPUTIME ((double) (clock ( )) / CLOCKS_PER_SEC)
#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#define PERM(j) (Lperm != NULL ? Lperm[j] : j)
#ifndef SPUMONI
#define SPUMONI 0 /* getting spumoni at run-time takes way too much time */
#endif
//#include "InvalidOperationException.h"
//#include "NullPointerException.h"

//#if !defined(_WIN32)
//#define dsymv dsymv_ /* Linux adds an underscore to command names */
//#define dgemv dgemv_
//#endif

/* This file pointer is used for the mread and mwrite mexFunctions.  It must
 * be a global variable, because the file pointer is not passed to the
 * sputil_error_handler function when an error occurs. */
//FILE *sputil_file = NULL ;

class RecoverMoments {
    
    //friend class PwgOptimiser;
    
private:
    
    Eigen::SparseMatrix<double> Y;
    Eigen::VectorXd y;
    
    cholmod_factor *compute_cholesky(
        cholmod_sparse *A,
        cholmod_common *cm);
        
    cholmod_sparse *compute_spinv(
            cholmod_factor *L,
            cholmod_common *cm);
    
    cholmod_sparse *cholmod_spinv(
            cholmod_factor *L,
            cholmod_common *cm);
    
    cholmod_sparse *cholmod_spinv_simplicial(
            cholmod_factor *L,
            cholmod_common *cm);
    
    cholmod_sparse *cholmod_spinv_super(
        cholmod_factor *L,
        cholmod_common *cm);
    
    cholmod_sparse *EigenDenseToCholmodSparseCopy(
            const Eigen::MatrixXd &in,
            cholmod_common *cm,
            double eps);
    
    cholmod_dense *EigenDenseToCholmodDenseCopy(
            const Eigen::MatrixXd &in,
            cholmod_common *cm);
    
    Eigen::MatrixXd CholmodDenseToEigenDenseCopy(
            cholmod_dense *in,
            cholmod_common *cm);
    
    cholmod_sparse *EigenSparseToCholmodSparseCopy(
            const Eigen::SparseMatrix<double> &in,
            cholmod_common *cm,
            double eps);
    
    Eigen::SparseMatrix<double> CholmodSparseToEigenSparseCopy(
            const cholmod_sparse *in,
            cholmod_common *cm,
            double eps);
    
    cholmod_dense *solve_system (
            cholmod_sparse *A,
            cholmod_factor *L,
            cholmod_dense *b,
            cholmod_common *cm) ;
    
public:
    
    //RecoverMoments();
    RecoverMoments( const Eigen::SparseMatrix<double> &A,
            const Eigen::VectorXd &b);
    ~RecoverMoments();
    Eigen::SparseMatrix<double> P;
    Eigen::VectorXd x;
    
};

extern "C" {

    //cholmod_sparse *cholmod_l_spinv( cholmod_factor *L, cholmod_common *Common ) ;
    
    void sputil_config (
            Int spumoni,
            cholmod_common *cm);
    
    void sputil_error_handler (
            int status,
            const char *file,
            int line,
            const char *message);
    
    //void spinv(
     //       Eigen::SparseMatrix<double> &P,
     //       Eigen::SparseMatrix<double> &Y);
    
    cholmod_sparse* sputil_get_sparse(
            const mxArray *Amatlab, /* MATLAB version of the matrix */
            cholmod_sparse *A,	    /* CHOLMOD version of the matrix */
            double *dummy,	    /* a pointer to a valid scalar double */
            Int stype);		    /* -1: lower, 0: unsymmetric, 1: upper */
    
    void cumsum2 (
            mwIndex *p,
            mwIndex *c,
            mwSize n);
};

mxArray *eigen_get_matlab_sparse(Eigen::SparseMatrix<double> &Y);
void spprint(const mxArray *mx);

#endif /* RECOVERMOMENTS_H */

/* 
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */
