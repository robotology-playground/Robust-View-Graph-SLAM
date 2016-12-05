/**
 * @file xxxx.cpp
 * @brief .
 * @detail .
 * @copyright Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
 * @authors Tariq Abuhashim, Nicolo' Genesio
 * @email t.abuhashim@gmail.com, nicogene@hotmail.it
 * @date Nov 2016
 * @acknowledgement This research has received funding from the European Union’s 
 * Seventh Framework Programme for research, technological development and demonstration 
 * under grant agreement No. 611909(KoroiBot).
 * @license Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
 */

#include "RecoverMoments.h"

// /*=====================*/
// RecoverMoments::RecoverMoments(void){
//     //std::cout << "Object is being created" << std::endl;
// }

/*=====================*/
RecoverMoments::RecoverMoments(const Eigen::SparseMatrix<double> &Y,
        const Eigen::VectorXd &y){
    
    /* initialise cholmod */
    cholmod_common Common, *cm ;
    cm = &Common;
    cholmod_l_start (cm) ;
    sputil_config (SPUMONI, cm) ;
    
    /* We have to find first the Cholesky factorization. */
    cholmod_sparse *A ;
    cholmod_factor *L ;
    A = EigenSparseToCholmodSparseCopy (Y, cm, 1e-12) ;
    L = compute_cholesky (A, cm) ;
    //printf("L->minor = %d\n", (int) L->minor);
    
    P.resize(y.size(),y.size()) ;
    x.resize(y.size(),1) ;
    x.setZero(y.size(),1) ;
    
    if ( L == NULL ) /* not positive definite - the expensive option */
    {
        std::cerr << "not positive definite - an alternative option" << std::endl ;
    }
    else /* positive definite - the cheap option */
    {
        /* Solve for x */
        cholmod_dense *b, *X;
        //b = cholmod_l_ones (A->nrow, 1, A->xtype, cm) ; // b = ones(n,1)
        b = EigenDenseToCholmodDenseCopy (y, cm) ;
        X = solve_system (A, L, b, cm) ;
        x = CholmodDenseToEigenDenseCopy (X, cm) ;
        
        /* Solve for P */
        cholmod_sparse *invA;
        invA = compute_spinv (L, cm) ;
        std::cout << "invA-computed" << std::endl ;
        sleep(10) ;
        //invA = cholmod_spinv (L, cm) ;
        P = CholmodSparseToEigenSparseCopy (invA, cm, 1e-12) ;
        std::cout << "P-computed" << std::endl ;
        //sleep(10) ;
        
        /* free workspace */
        cholmod_l_free_sparse (&invA, cm) ;
        cholmod_l_free_dense (&X, cm) ;
        cholmod_l_free_dense (&b, cm) ;
        std::cout << "memory 1 removed" << std::endl ;
        //sleep(10) ;
    }
    
    /* free workspace and the CHOLMOD L */
    cholmod_l_free_factor (&L, cm) ;                    /* free matrices */
    cholmod_l_free_sparse (&A, cm) ;
    cholmod_l_finish (cm) ;                            /* finish CHOLMOD */
    std::cout << "memory 2 removed" << std::endl ;
    //sleep(10) ;
    //delete cm ;
}

/*=====================*/
RecoverMoments::~RecoverMoments(void){
    //std::cout << "RecoverMoments object is being deleted" << std::endl;
}

/*=====================*/
/* Define function pointers and other parameters for a mexFunction */
void sputil_config (Int spumoni, cholmod_common *cm)
{
    
    /* convert to packed LDL' when done */
    cm->final_asis = FALSE ;
    cm->final_super = FALSE ;
    cm->final_ll = FALSE ;
    cm->final_pack = TRUE ;
    cm->final_monotonic = TRUE ;
    
    /* since numerically zero entries are NOT dropped from the symbolic
     * pattern, we DO need to drop entries that result from supernodal
     * amalgamation. */
    cm->final_resymbol = TRUE ;
    cm->quick_return_if_not_posdef = TRUE ;
    
    /* cholmod_l_solve must return a real or zomplex X for MATLAB */
    cm->prefer_zomplex = TRUE ;
    
    //cm->print = 5 ;
    
    //cm->supernodal = CHOLMOD_SIMPLICIAL ;
    //cm->supernodal = CHOLMOD_SUPERNODAL ;
    
    /* use mxMalloc and related memory management routines */
    //cm->malloc_memory  = mxMalloc ;
    //cm->free_memory    = mxFree ;
    //cm->realloc_memory = mxRealloc ;
    //cm->calloc_memory  = mxCalloc ;
    
    /* printing and error handling */
    if (spumoni == 0)
    {
        /* do not print anything from within CHOLMOD */
        cm->print = -1 ;
        //cm->print_function = NULL ;
    }
    else
    {
        /* spumoni = 1: print warning and error messages.  cholmod_l_print_*
         *	routines will print a one-line summary of each object printed.
         * spumoni = 2: also print a short summary of each object.
         */
        cm->print = spumoni + 2 ;
        //cm->print_function = mexPrintf ;
    }
    
    /* error handler */
    //cm->error_handler  = sputil_error_handler ;
    
    /* complex arithmetic */
    //cm->complex_divide = cholmod_l_divcomplex ;
    //cm->hypotenuse     = cholmod_l_hypot ;
    
#ifndef NPARTITION
#if defined(METIS_VERSION)
#if (METIS_VERSION >= METIS_VER(4,0,2))
    /* METIS 4.0.2 uses function pointers for malloc and free */
    METIS_malloc = cm->malloc_memory ;
    METIS_free   = cm->free_memory ;
#endif
#endif
#endif
    
    /* Turn off METIS memory guard.  It is not needed, because mxMalloc will
     * safely terminate the mexFunction and free any workspace without killing
     * all of MATLAB.  This assumes cholmod_make was used to compile CHOLMOD
     * for MATLAB. */
    //cm->metis_memory = 0.0 ;
}

/*=====================*/
cholmod_factor *RecoverMoments::compute_cholesky(
        cholmod_sparse *A,
        cholmod_common *cm)
{
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "compute_cholesky: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        return NULL ;
    }
    
    double t;
    t = CPUTIME ;
    
    if (A == NULL || A->stype == 0)               /* A must be symmetric */
    {
        return NULL ;
    }
    
    /* analyze and factorise */
    cholmod_factor *L;
    L = cholmod_l_analyze (A, cm) ;
    if (A->stype == 0)
    {   //printf ("Factorising A*A'+beta*I\n") ;
        double beta [2] ;
        beta [0] = 0 ;
        beta [1] = 0 ;
        cholmod_l_factorize_p (A, beta, NULL, 0, L, cm) ;
    }
    else
    {   //printf ("Factorising A\n") ;
        cholmod_l_factorize (A, L, cm) ;
    }
    if (cm->status != CHOLMOD_OK)
    {
        std::cerr << "matrix is not positive definite" << std::endl ;
        return NULL ;
    }
    
    t = CPUTIME - t ;
    t = MAX (t, 0) ;
    //std::cout << "compute_cholesky: " << t << std::endl ;
    
    return L ;
}

/*=====================*/
cholmod_sparse* RecoverMoments::compute_spinv(cholmod_factor *L,
        cholmod_common *cm)
{
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "compute_spinv: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return NULL ;
    }
    
    double t ;
    t = CPUTIME ;
    
    int sorted = true, packed = true, stype = -1 ; // NEED TO MAKE THE MATRIX SYMMETRIC
    
    /* ----------------------------------------------------------------- */
    /* convert L to a sparse matrix */
    /* ----------------------------------------------------------------- */
    
    cholmod_sparse *Lsparse ;
    Lsparse = cholmod_l_factor_to_sparse (L, cm) ;
    if (Lsparse->xtype == CHOLMOD_COMPLEX)
    {
        //mexErrMsgTxt ("matrix is complex") ;
    }
    mwSize nnz = (mwSize)(Lsparse->nzmax);
    mwIndex *I, *J;
    I = (mwIndex *)(Lsparse->i) ; // A->i, an integer array of size A->nzmax.
    J = (mwIndex *)(Lsparse->p) ; // A->p, an integer array of size A->ncol+1.
    double *C = (double*)(Lsparse->x) ; // A->x, a double array of size A->nzmax or twice that for the complex case.
    size_t ncol = (size_t)L->n;
    
    std::cout << "spinv step #1" << std::endl ;
    //sleep(10) ;
    
    /* ----------------------------------------------------------------- */
    /* Evaluate the sparse inverse */
    /* ----------------------------------------------------------------- */
    
    //char *uplo = "L", *trans = "N";
    double done = 1.0, dzero = 0.0;
    const ptrdiff_t one = 1;
    mwIndex l, k2, h, k, i, j, ik, lfi;
    std::ptrdiff_t lfi_si;

    //double *fil=(double*)(calloc((mwSize)1,sizeof(double)));
    //double *zt=(double*)(mxCalloc((mwSize)1,sizeof(double)));
    //double *Zt=(double*)(mxCalloc((mwSize)1,sizeof(double)));
    //double *zz=(double*)(mxCalloc((mwSize)1,sizeof(double)));     
    double *fil = new double [(int)1]();
    double *zt = new double [(int)1]();
    double *Zt = new double [(int)1]();
    double *zz = new double [(int)1]();

    C[nnz-1] = 1.0/C[J[ncol-1]];/*set the last element of sparse inverse*/
    for (int j=ncol-2; j!=-1; j--)
    {
        lfi = J[j+1]-(J[j]+1);
        
        /* if (lfi > 0) */
        if (J[j+1]-(J[j]+1) > 0)
        {
            /*	printf("lfi = %u \n ", lfi);
             * printf("lfi*double = %u \n", (mwSize)lfi*sizeof(double));
             * printf("lfi*lfi*double = %u \n", (mwSize)lfi*(mwSize)lfi*sizeof(double));
             * printf("\n \n");
             */
            //std::cout << lfi << std::endl;
            delete[] fil ;
            delete[] zt ;
            delete[] Zt ;
            fil = new double [lfi](); /* memory for the sparse inverse elements to be evaluated */
            zt = new double [lfi](); /* memory for the needed sparse inverse elements */
            Zt = new double [lfi*lfi]();
            
            /* Set the lower triangular for Zt */
            /* take the j'th lower triangular column of the Cholesky */
            for (i=0; i<lfi; i++) fil[i] = C[J[j]+i+1];
            k2 = 0 ; 
            for (k=J[j]+1; k<J[j+1]; k++){
                ik = I[k];
                h = k2;
                for (l=J[ik]; l<=J[ik+1]; l++){
                    //if (I[l] == I[ J[j]+h+1 ]){
                    //    Zt[h+lfi*k2] = C[l];
                    //    h++;
                    //}
                }
                k2++;
            }
            /* evaluate zt = fil*Zt */
            lfi_si = (mwSignedIndex) lfi;
            dsymv_("L", &lfi_si, &done, Zt, &lfi_si, fil, &one, &dzero, zt, &one);
            /* Set the evaluated sparse inverse elements, zt, into C */
            k=lfi-1;
            for (i = J[j+1]-1; i!=J[j] ; i--){
                C[i] = -zt[k];
                k--;
            }
            /* evaluate the j'th diagonal of sparse inverse */
            dgemv_("N", &one, &lfi_si, &done, fil, &one, zt, &one, &dzero, zz, &one);
            C[J[j]] = 1.0/C[J[j]] + zz[0];
       }
       else
       {
           /* evaluate the j'th diagonal of sparse inverse */
          C[J[j]] = 1.0/C[J[j]];
       }
        
    }
    /* Free the temporary variables */
    delete[] fil ; fil = NULL ;
    delete[] zt; zt = NULL ;
    delete[] Zt; Zt = NULL ;
    delete[] zz; zz = NULL ;
        
    std::cout << "spinv step #2" << std::endl ;
    //sleep(10) ;
    
    /* ----------------------------------------------------------------- */
    /* Permute the elements according to r(q) = 1:n */
    /* Done only if the Cholesky was evaluated here */
    /* ----------------------------------------------------------------- */
    
    cholmod_sparse* Bm;
    Bm = cholmod_l_allocate_sparse ( ncol, ncol, nnz, sorted, packed, stype,
            CHOLMOD_REAL, cm ) ;
    //nnz = Lsparse->nzmax ; // Number of non-zero elements
    //int xtype = L->xtype ; // Result has the same xtype as the factor
    //Bm = cholmod_spzeros (ncol, ncol, nnz, CHOLMOD_REAL, cm) ;
    
    if (Bm == NULL){
        //mexErrMsgTxt ("memory #1 failed") ;
        std::cerr << "memory #1 failed" << std::endl;
        //cholmod_l_free_sparse (&Bm, cm) ;
        //cholmod_l_free_sparse (&Lsparse, cm) ;
        return 0;
        //throw InvalidOperationException("cholmod_l_allocate_sparse failed",
        //    __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    mwIndex *It = (mwIndex*)(Bm->i);
    mwIndex *Jt = (mwIndex*)(Bm->p);
    double *Ct = (double*)(Bm->x); /*Ct = C(r,r)*/
    
    mwIndex *r, *w, *w2;
    r = (mwIndex *)(L->Perm) ; /* fill reducing ordering */
    /* mxCalloc initializes this newly allocated memory to 0 */
    /* new mwIndex [ncol]() : dynamically allocates space for ncol elements 
     * of type mwIndex, () initialises with 0 */
    w = new mwIndex [ncol]() ;  
    
    /* column counts of Am */
    /* count entries in each column of Bm */
    for (j=0; j<ncol; j++){
        k = r ? r[j] : j ;           /* column j of Bm is column k of Am */
        for (l=J[j] ; l<J[j+1] ; l++){
            i = I[l];
            ik = r ? r[i] : i ;           /* row i of Bm is row ik of Am */
            w[ MAX(ik,k) ]++;
        }
    }
    cumsum2 (Jt, w, ncol) ;
    for (j=0; j<ncol; j++){
        k = r ? r[j] : j ;           /* column j of Bm is column k of Am */
        for (l=J[j] ; l<J[j+1] ; l++){
            i= I[l];
            ik = r ? r[i] : i ;           /* row i of Bm is row ik of Am */
            It [k2 = w[MAX(ik,k)]++ ] = MIN(ik,k);
            Ct[k2] = C[l];
        }
    }
    delete[] w ;
        
    std::cout << "spinv step #3" << std::endl ;
    //sleep(10) ;
    
    /* ----------------------------------------------------------------- */
    /* Transpose the permuted (upper triangular) matrix
     * Bm into Am (this way we get sorted columns) */
    /* ----------------------------------------------------------------- */
    
    w = new mwIndex [ncol] () ;
    for (i=0 ; i<Jt[ncol] ; i++) w[It[i]]++;         /* row counts of Bm */
    cumsum2 (J, w, ncol) ;                               /* row pointers */
    for (j=0 ; j<ncol ; j++){
        for (i=Jt[j] ; i<Jt[j+1] ; i++){
            I[ l=w[ It[i] ]++ ] = j;
            C[l] = Ct[i];
        }
    }
    delete[] w ;
    cholmod_l_free_sparse (&Bm, cm) ;
        
    std::cout << "spinv step #4" << std::endl ;
    //sleep(10) ;
    
    /* ----------------------------------------------------------------- */
    /* Fill the upper triangle of the sparse inverse */
    /* ----------------------------------------------------------------- */

   w = new mwIndex [ncol]() ; 
   w2 = new mwIndex [ncol]() ; 
    for (k=0; k<J[ncol]; k++) 
        w[I[k]]++;/* row counts of the lower triangular */
    for (k=0; k<ncol; k++) 
        w2[k] = w[k] + J[k+1] - J[k] - 1;
    /* column counts of the sparse inverse */
    nnz = (mwSize)2*nnz - ncol;              /* The number of nonzeros in Z */
    Bm = cholmod_l_allocate_sparse ( ncol, ncol, nnz, sorted, packed, stype,
            CHOLMOD_REAL, cm ) ;
    //Bm = cholmod_spzeros (ncol, ncol, nnz, CHOLMOD_REAL, cm) ;
    //cholmod_dense *A;
    //A = cholmod_zeros (ncol, ncol, CHOLMOD_REAL, cm) ;
    //Bm = cholmod_dense_to_sparse(A, TRUE, cm) ;
    //Bm->stype = 1 ; // NEED TO MAKE THE MATRIX SYMMETRIC
    
    if (Bm == NULL){
        std::cerr << "memory failed" << std::endl;
        //throw InvalidOperationException("cholmod_l_allocate_sparse failed",
        //    __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    It = (mwIndex*)(Bm->i);
    Jt = (mwIndex*)(Bm->p);
    Ct = (double*)(Bm->x);
    cumsum2(Jt, w2, ncol);                        /* column starting points */
    for (j = 0 ; j < ncol ; j++){              /* fill the upper triangular */
        for (k = J[j] ; k < J[j+1] ; k++){
            It[l = w2[ I[k]]++] = j ;	/* place C(i,j) as entry Ct(j,i) */
            if (Ct) Ct[l] = C[k] ;
        }
    }
    for (j = 0 ; j < ncol ; j++){              /* fill the lower triangular */
        for (k = J[j]+1 ; k < J[j+1] ; k++){
            It[l = w2[j]++] = I[k] ;    /* place C(j,i) as entry Ct(j,i) */
            if (Ct) Ct[l] = C[k] ;
        }
    }
    
    cholmod_l_free_sparse (&Lsparse, cm) ;
    delete[] w ; w = NULL ;
    delete[] w2 ; w2 = NULL ;

    t = CPUTIME - t ;
    t = MAX (t, 0) ;
    
    //std::cout << "compute_spinv: " << t << std::endl ;
    
    return Bm ;
        
    std::cout << "spinv step #5" << std::endl ;
    //sleep(10) ;
    
}

/*=====================*/
void cumsum2 (mwIndex *p, mwIndex *c, mwSize n){
    int i;
    int nz = 0;
    if(!p || !c) return;
    for (i=0;i<n;i++){
        p[i]=nz;
        nz+=c[i];
        c[i]=p[i];
    }
    p[n]=nz;
}

/*=====================*/
cholmod_sparse* RecoverMoments::EigenDenseToCholmodSparseCopy(const Eigen::MatrixXd& in,
        cholmod_common* cm, double eps) {
    
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "EigenDenseToCholmodSparseCopy: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return (NULL) ;
    }
    
    //if (cholmod == NULL)
    //    throw NullPointerException("cholmod", __FILE__, __LINE__,
    //            __PRETTY_FUNCTION__);
    size_t nzmax = 0;
    for (std::ptrdiff_t i = 0; i < in.rows(); ++i)
        for (std::ptrdiff_t j = 0; j < in.cols(); ++j)
            if (std::fabs(in(i, j)) > eps)
                nzmax++;
    int sorted = true, packed = true, stype = -1;
    cholmod_sparse* out = cholmod_l_allocate_sparse ( in.rows(), in.cols(),
            nzmax, sorted, packed, stype, CHOLMOD_REAL, cm ) ;
    if (out == NULL){
        std::cerr << "memory failed" << std::endl;
        //return 0;
        //throw InvalidOperationException("cholmod_l_allocate_sparse failed",
        //    __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    double* values = (double*)(out->x);
    std::ptrdiff_t* col_ptr = (std::ptrdiff_t*)(out->p);
    std::ptrdiff_t* row_ind = (std::ptrdiff_t*)(out->i);
    std::ptrdiff_t rowIt = 0, colIt = 1;
    for (std::ptrdiff_t c = 0; c < in.cols(); ++c) {
        for (std::ptrdiff_t r = 0; r < in.rows(); ++r)
            if (std::fabs(in(r, c)) > eps) {
                values[rowIt] = in(r, c);
                row_ind[rowIt] = r;
                rowIt++;
            }
        col_ptr[colIt] = rowIt;
        colIt++;
    }
    out->itype = CHOLMOD_LONG ;
    out->dtype = CHOLMOD_DOUBLE ;
    out->stype = -1 ;
    return out;
}

/*=====================*/
cholmod_dense* RecoverMoments::EigenDenseToCholmodDenseCopy(const Eigen::MatrixXd& in,
        cholmod_common* cm) {
    
    cholmod_dense *out = (NULL) ;
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "EigenDenseToCholmodDenseCopy: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return out ;
    }
    size_t nrow, ncol;
    nrow = in.rows();
    ncol = in.cols();
    out = cholmod_l_allocate_dense ( nrow, ncol, nrow*ncol, CHOLMOD_REAL, cm ) ;
    double *y;
    y = (double *) out->x;
    for(size_t j = 0; j < ncol; j++)// Column
        for(size_t i = 0; i < nrow; i++)// Row
            y[i+j*nrow] = in.coeffRef(i,j);
    
    return out;
}

/*=====================*/
Eigen::MatrixXd RecoverMoments::CholmodDenseToEigenDenseCopy(cholmod_dense* in,
        cholmod_common* cm) {
    
    size_t ncol, nrow ;
    ncol = (size_t) in->ncol ;
    nrow = (size_t) in->nrow ;
    Eigen::MatrixXd out (nrow,ncol) ;
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "CholmodDenseToEigenDenseCopy: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return out ;
    }
    
    double *y ;
    y = (double *) in->x ;
    for(size_t j = 0; j < ncol; j++)// Column
        for(size_t i = 0; i < nrow; i++)// Row
            out.coeffRef(i,j) = y[i+j*nrow] ;
    
    return out ;
}

/*=====================*/
cholmod_sparse* RecoverMoments::EigenSparseToCholmodSparseCopy(
        const Eigen::SparseMatrix<double>& in,
        cholmod_common* cm, double eps) {
    
    cholmod_sparse *out = (NULL) ;
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "EigenSparseToCholmodSparseCopy: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return out ;
    }
    double t;
    t = CPUTIME ;
    
    //if (cholmod == NULL)
    //    throw NullPointerException("cholmod", __FILE__, __LINE__,
    //            __PRETTY_FUNCTION__);
    size_t nzmax = 0;
    for (std::ptrdiff_t i = 0; i < in.outerSize(); ++i)
        for (Eigen::SparseMatrix<double>::InnerIterator it(in,i); it; ++it)
            if (std::fabs(it.value()) > eps)
                nzmax++;
    int sorted = true, packed = true, stype = -1;
    out = cholmod_l_allocate_sparse ( in.outerSize(), in.outerSize(),
            nzmax, sorted, packed, stype, CHOLMOD_REAL, cm ) ;
    if (out == NULL){
        std::cerr << "memory failed" << std::endl;
        return 0;
        //throw InvalidOperationException("cholmod_l_allocate_sparse failed",
        //    __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    double* values = (double*)(out->x);
    std::ptrdiff_t* col_ptr = (std::ptrdiff_t*)(out->p);
    std::ptrdiff_t* row_ind = (std::ptrdiff_t*)(out->i);
    std::ptrdiff_t rowIt = 0, colIt = 1;
    for (std::ptrdiff_t i = 0; i<in.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(in,i); it; ++it)
            if (std::fabs(it.value()) > eps) {
                values[rowIt] = it.value();
                row_ind[rowIt] = it.row();
                rowIt++;
            }
        col_ptr[colIt] = rowIt;
        colIt++;
    }
    out->itype = CHOLMOD_LONG ;
    out->dtype = CHOLMOD_DOUBLE ;
    out->stype = -1 ; /* use lower part of A */
    //out->itype = CHOLMOD_INT;
    //cholmod_print_sparse (out, "A", cm);
    
    t = CPUTIME - t ;
    t = MAX (t, 0) ;
    //printf ("EigenSparseToCholmodSparseCopy:  %12.4f\n", t) ;
    
    return out;
}

/*=====================*/
Eigen::SparseMatrix<double> RecoverMoments::CholmodSparseToEigenSparseCopy(
        const cholmod_sparse* in, cholmod_common* cm, double eps) {
    
    Eigen::SparseMatrix<double> out (in->nrow, in->ncol) ;
    
    //std::cout << "dims " << in->nrow << " " << in->ncol << std::endl ;
    //std::cout << "status " << cm->status << std::endl ;
    
    if (cm->status < CHOLMOD_OK)
    {
        std::cout << "CholmodSparseToEigenSparseCopy: " ;
        std::cerr << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
        //return out ;
    }
    double* C = (double*)(in->x);
    std::ptrdiff_t* J = (std::ptrdiff_t*)(in->p);
    std::ptrdiff_t* I = (std::ptrdiff_t*)(in->i);
    
    std::vector<T> tripletList;
    for (int j=0; j<(in->ncol); j++)             /* Loop through columns */
        for (int i=J[j]; i<J[j+1]; i++) /* Loop through non-zeros in ith column */
            if (std::fabs(C[i])>eps) tripletList.push_back(T(I[i], j, C[i]));
    out.setFromTriplets(tripletList.begin(), tripletList.end());
    
    //std::cout << out << std::endl ;
    
    //delete[] C ; // causes error when emptying invA
    //delete[] I ; // causes error when emptying invA
    //delete[] J ; // causes error when emptying invA
    return out;
}

// /*=====================*/
// mxArray *eigen_get_matlab_sparse(
//         Eigen::SparseMatrix<double> &Y){
//     /* MATLAB Sparse to Cholmod Sparse */
//     /* HOW TO USE ? - EXAMPLE */
//     //mxArray *Amatlab;
//     //cholmod_sparse Amatrix ;
//     //Amatlab = eigen_get_matlab_sparse (Y);
//     //double dummy = 0;
//     //A = sputil_get_sparse (Amatlab, &Amatrix, &dummy, -1) ;
//     
//     /* get sparse matrix A, use tril(A)  */
//     //cholmod_sparse *A;	    /* CHOLMOD version of the matrix */
//     long unsigned int nzmax = Y.nonZeros();
//     mxArray *Amatlab = mxCreateSparse(Y.outerSize(), Y.outerSize(), nzmax, mxREAL);
//     double *pr = mxGetPr(Amatlab); /* Non-zero elements */
//     long unsigned int *ir = mxGetIr(Amatlab); /* Row indexing */
//     long unsigned int *jc = mxGetJc(Amatlab); /* Column count */
//     int index=0, kk;
//     for (kk=0; kk < Y.outerSize(); ++kk)
//     {
//         jc[kk] = index;
//         for (Eigen::SparseMatrix<double>::InnerIterator it(Y,kk); it; ++it)
//         {
//             ir[index] = it.row();
//             pr[index] = it.value();
//             index++;
//         }
//     }
//     jc[Y.outerSize()] = index;
//     //spprint(Amatlab); /* prints out the sparse MATLAB matrix */
//     return (Amatlab) ;
// }

// /*=====================*/
// void spprint(const mxArray *mx)
// {
//     
//     mwSize n, nrow;
//     mwIndex *ir, *jc;
//     mwIndex j, x, y;
//     double *pr;
//     if( !mxIsSparse(mx) ) return;
//     n = mxGetN(mx);
//     pr = mxGetPr(mx);
//     ir = mxGetIr(mx);
//     jc = mxGetJc(mx);
//     for( y=0; y<n; y++ ) {
//         nrow = jc[y+1] - jc[y];
//         for( x=0; x<nrow; x++ ) {
//             mexPrintf("   (%d,%d)    %g\n",(*ir++)+1,y+1,*pr++);
//         }
//     }
// }

// /*=====================*/
// /* returns the sparse solution */
// cholmod_sparse *RecoverMoments::cholmod_spinv(
//         cholmod_factor *L,	/* factorization to use */
//         cholmod_common *cm)
// {
//     
//     if (cm->status < CHOLMOD_OK)
//     {
//         std::cout << "cholmod_spinv: " ;
//         std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
//         return (NULL) ;
//     }
//     
//     double t;
//     t = CPUTIME ;
//     cholmod_sparse* Bm;
//     //ASSERT (L->xtype != CHOLMOD_PATTERN) ;  /* L is not symbolic */
//     
//     /* ---------------------------------------------------------------------- */
//     /* check inputs */
//     /* ---------------------------------------------------------------------- */
//     
//     //RETURN_IF_NULL_COMMON (NULL) ;
//     //RETURN_IF_NULL (L, NULL) ;
//     //RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, NULL) ;
//     //cm->status = CHOLMOD_OK ;
//     
//     /*
//      * Compute the sparse inverse.
//      */
//     if (L->is_super)
//     {
//         std::cout << "spinv_super" << std::endl;
//         Bm = cholmod_spinv_super (L, cm) ;
//         t = CPUTIME - t ;
//         t = MAX (t, 0) ;
//         printf ("cholmod_spinv_super:  %12.4f\n", t) ;
//         return 0;
//     }
//     else
//     {
//         std::cout << "spinv_simplicial" << std::endl;
//         Bm = cholmod_spinv_simplicial (L, cm) ;
//         t = CPUTIME - t ;
//         t = MAX (t, 0) ;
//         printf ("cholmod_spinv_simplicial:  %12.4f\n", t) ;
//         return Bm;
//     }
// }

// /*=====================*/
// /* returns the sparse solution */
// cholmod_sparse *RecoverMoments::cholmod_spinv_super(
//         cholmod_factor *L,	/* factorization to use */
//         cholmod_common *cm)
// {
//     if (cm->status < CHOLMOD_OK)
//     {
//         std::cout << "cholmod_spinv_super: " ;
//         std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
//         return NULL ;
//     }
//     
// }

// /*=====================*/
// /* returns the sparse solution */
// cholmod_sparse *RecoverMoments::cholmod_spinv_simplicial(
//         cholmod_factor *L,	/* factorization to use */
//         cholmod_common *cm)
// {
//     if (cm->status < CHOLMOD_OK)
//     {
//         std::cout << "cholmod_spinv_simplicial: " ;
//         std::cout << "NULL Common, out of memory, or inputs invalid" << std::endl ;
//         return (NULL) ;
//     }
//     
//     int xtype ;
//     cholmod_sparse *X ;
//     double *Lx, *Lz, *Xx, *Xz, *V, *z, *Lxj ;
//     double djj ;
//     Int *Li, *Lp, *Xp, *Xi ;
//     Int *perm, *Lperm, *ncol ;
//     Int n, kmin, kmax, iz, jz, il, jl, kl, ix, jx, kx, ip, jp ;
//     size_t nz, maxsize ;
//     
//     mwSignedIndex nj;
//     double minus_one[2], zero[2] ;
//     //double one[2] ;
//     minus_one[0] = -1.0 ;
//     minus_one[1] = 0.0 ;
//     //one[0] = 1.0 ;
//     //one[1] = 0.0 ;
//     zero[0] = 0.0 ;
//     zero[1] = 0.0 ;
//     
//     double done = 1.0, dzero = 0.0;
//     const ptrdiff_t one = 1;
//     
//     n = L->n ; // Dimensionality of the matrix
//     nz = L->nzmax ; // Number of non-zero elements
//     xtype = L->xtype ; // Result has the same xtype as the factor
//     
//     // Allocate the result X
//     X = cholmod_spzeros (n, n, nz, xtype, cm) ;
//     if (cm->status < CHOLMOD_OK)
//     {
//         std::cout << "memory #3 failed" << std::endl;
//         cholmod_free_sparse (&X, cm) ;
//         return NULL ;
//     }
//     
//     // Shorthand notation
//     Xp = (long int*)(X->p) ;
//     Xi = (long int*)(X->i) ;
//     Xx = (double*)(X->x) ;
//     Xz = (double*)(X->z) ;
//     Lp = (long int*)(L->p) ;
//     Li = (long int*)(L->i) ;
//     Lx = (double*)(L->x) ;
//     Lz = (double*)(L->z) ;
//     Lperm = (long int*)(L->Perm) ;
//     
//     V = NULL ;
//     z = NULL ;
//     Lxj = NULL ;
//     
//     /*
//      * Compute the mapping to the permuted result:
//      * X->x[perm[i]] ~ L->x[i]
//      * Both X and L are lower triangular
//      */
//     
//     /* Count non-zeros on columns */
//     for (jl = 0; jl < n; jl++)
//     {
//         jp = PERM(jl) ; // permuted column
//         for (kl = Lp[jl]; kl < Lp[jl+1]; kl++)
//         {
//             il = Li[kl] ;    // row of L
//             ip = PERM(il) ; // permuted row
//             // Increase the number of elements in the column
//             jx = MIN(ip,jp) ;
//             Xp[jx+1]++ ;
//         }
//     }
//     
//     /* Compute column pointers by computing cumulative sum */
//     maxsize = 0 ;
//     for (jx = 1; jx <= n; jx++)
//     {
//         if (Xp[jx] > maxsize)
//             maxsize = Xp[jx] - 1 ; // number of non-zeros (without diagonal)
//         Xp[jx] += Xp[jx-1] ;
//     }
//     
//     /* Add row indices */
//     ncol = (Int*)calloc(n, sizeof(Int)) ;
//     perm = (long int*)malloc(nz*sizeof(Int)) ; // permutation mapping
//     
//     for (jl = 0; jl < n; jl++)
//     {
//         jp = PERM(jl) ; // permuted column
//         for (kl = Lp[jl]; kl < Lp[jl+1]; kl++)
//         {
//             il = Li[kl] ;          // row of L
//             ip = PERM(il) ;       // permuted row
//             jx = MIN(ip,jp) ;      // column of X
//             ix = MAX(ip,jp) ;      // row of X
//             kx = Xp[jx]+ncol[jx] ; // index of X
//             // Increase elements on the column
//             ncol[jx]++ ;
//             // Add row index
//             Xi[kx] = ix ;
//             // Mapping X[perm[i]] ~ L[i]
//             perm[kl] = kx ;
//         }
//     }
//     free(ncol) ;
//     X->sorted = FALSE ;
//     
//     // Allocate memory for a temporary matrix and vector
//     z = (double*)malloc((maxsize+1)*sizeof(double)) ;
//     V = (double*)malloc((maxsize*maxsize)*sizeof(double)) ;
//     
//     if (L->is_ll)
//     {
//         
//         switch (xtype)
//         {
//             case CHOLMOD_REAL:
//                 //ERROR (CHOLMOD_INVALID,"Real xtype for L*L' not implemented.") ;
//                 break ;
//                 
//             case CHOLMOD_COMPLEX:
//                 //ERROR (CHOLMOD_INVALID,"Complex xtype for L*L' not implemented.") ;
//                 break ;
//                 
//             case CHOLMOD_ZOMPLEX:
//                 //ERROR (CHOLMOD_INVALID,"Zomplex xtype for L*L' not implemented.") ;
//                 break ;
//                 
//         }
//     }
//     else
//     {
//         switch (xtype)
//         {
//             case CHOLMOD_REAL:
//                 
//                 for (jl = n-1; jl >= 0; jl--)
//                 {
//                     // Indices of non-zero elements in j-th column
//                     kmin = Lp[jl];         // first index
//                     kmax = Lp[jl+1] - 1;   // last index
//                     nj = kmax - kmin; // number of non-zero elements (without diagonal)
//                     
//                     // Diagonal entry of D: D[j,j]
//                     djj = Lx[kmin] ;
//                     if (kmax > kmin)
//                     {
//                         // j-th column vector of L (without the
//                         // diagonal element and zeros)
//                         Lxj = Lx + (kmin+1) ;
//                         
//                         // Form Z
//                         for (jz = 0; jz < nj; jz++)
//                         {
//                             // Row index of the (jz+1):th non-zero
//                             // element on column j
//                             // = relevant column index of X
//                             jx = Li[kmin+1+jz] ;
//                             // Index of the diagonal element on column
//                             // jx
//                             kx = Lp[jx] ;
//                             
//                             // Set lower triangular elements of column
//                             // jz (no need to set upper triangular
//                             // elements because of the symmetry).
//                             for (iz = jz; iz < nj; iz++)
//                             {
//                                 ix = Li[kmin+1+iz] ;
//                                 // Find X[row,jx]
//                                 while (Li[kx] < ix) kx++ ;
//                                 // Set Z[iz,jz] = X[ix,jx]
//                                 V[iz+jz*nj] = Xx[perm[kx]] ;
//                             }
//                         }
//                         
//                         // DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
//                         //BLAS_dsymv("L", nj, one, V, nj, Lxj, 1, zero, z, 1) ;
//                         dsymv_("L", &nj, &done, V, &nj, Lxj, &one, &dzero, z, &one) ;
//                         //dsymv_("L", &lfi_si, &done, Zt, &lfi_si, fil, &one, &dzero, zt, &one) ;
//                         
//                         // Copy the result to the lower part of X
//                         for (iz = 0; iz < nj; iz++)
//                         {
//                             kx = kmin + 1 + iz ;
//                             Xx[perm[kx]] = -z[iz] ;
//                         }
//                         
//                         // Compute the diagonal element X[j,j]
//                         // DDOT(N,DX,INCX,DY,INCY)
//                         //BLAS_ddot(nj, z, 1, Lxj, 1, Xx[perm[kmin]]) ;
//                         Xx[perm[kmin]] = ddot_(&nj, z, &one, Lxj, &one) ; // conflicted with matlab extras.
//                         Xx[perm[kmin]] += 1.0/djj ;
//                         
//                     }
//                     else
//                     {
//                         // Compute the diagonal element X[j,j]
//                         Xx[perm[kmin]] = 1.0/djj ;
//                     }
//                 }
//                 break ;
//                 
//             case CHOLMOD_COMPLEX:
//                 //ERROR (CHOLMOD_INVALID,"Complex xtype for L*D*L' not implemented.") ;
//                 break ;
//                 
//             case CHOLMOD_ZOMPLEX:
//                 //ERROR (CHOLMOD_INVALID,"Zomplex xtype for L*D*L' not implemented.") ;
//                 break ;
//         }
//         
//     }
//     
//     free(V) ;
//     free(z) ;
//     free(perm) ;
//     
//     if (cm->status == CHOLMOD_OK)
//     {
//         // The result is symmetric but only the lower triangular part was
//         // computed.
//         X->stype = -1 ;
//         
//         // Sort columns (is it necessary?)
//         cholmod_sort (X, cm) ;
//         
//         if (cm->status == CHOLMOD_OK)
//             return (X) ;
//     }
//     
//     cholmod_l_free_sparse (&X, cm) ;
//     return NULL ;
//     
// }

/*=====================*/
/* Solve for x */
cholmod_dense *RecoverMoments::solve_system ( cholmod_sparse *A,
        cholmod_factor *L, cholmod_dense *b, cholmod_common *cm)
{
    cholmod_dense *X, *r, *w;
    double one [2], zero [2], minusone [2], beta [2] ;
    zero [0] = 0 ; zero [1] = 0 ;
    one [0] = 1 ; one [1] = 0 ;
    minusone [0] = -1 ; minusone [1] = 0 ;
    beta [0] = 1e-6 ; beta [1] = 0 ;
    int i, n, isize, xsize, ordering, xtype ;
    xtype = A->xtype ;
    n = A->nrow ;
    double *Bx, *Rx, *Xx ;

    X = cholmod_l_solve (CHOLMOD_A, L, b, cm) ; // solves Ax=b
    
    Bx = (double*) b->x ;
    if (A->stype == 0)
    {
        /* (AA'+beta*I)x=b is the linear system that was solved */
        /* w = A'*X */
        w = cholmod_l_allocate_dense (A->ncol, 1, A->ncol, xtype, cm) ;
        cholmod_l_sdmult (A, 2, one, zero, X, w, cm) ;
        /* R = B - beta*X */
        r = cholmod_l_zeros (n, 1, xtype, cm) ;
        Rx = (double*) r->x ;
        Xx = (double*) X->x ;
        if (xtype == CHOLMOD_REAL)
        {
            for (i = 0 ; i < n ; i++)
            {
                Rx [i] = Bx [i] - beta [0] * Xx [i] ;
            }
        }
        else
        {
            /* complex case */
            for (i = 0 ; i < n ; i++)
            {
                Rx [2*i  ] = Bx [2*i  ] - beta [0] * Xx [2*i  ] ;
                Rx [2*i+1] = Bx [2*i+1] - beta [0] * Xx [2*i+1] ;
            }
        }
        /* R = A*W - R */
        cholmod_l_sdmult (A, 0, one, minusone, w, r, cm) ;
        cholmod_free_dense (&w, cm) ;
    }
    else
    {
        /* Ax=b was factorized and solved, R = B-A*X */
        r = cholmod_l_copy_dense (b, cm) ;
        cholmod_l_sdmult (A, 0, minusone, one, X, r, cm) ;
    }
    
    cholmod_l_free_dense (&r, cm) ;
    
    return X ;
}

/* 
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */
