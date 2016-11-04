#include "mex.h"
#include <iostream> // std::cout
#include <string>
#include <math.h>       /* floor */
#include <list>
#include <vector>
#include <memory.h>

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* create a pointer to the real data in the input matrix  */
    double *p = mxGetPr(prhs[0]); /* image points */
    double *dims = mxGetPr(prhs[1]);
    double *N = mxGetPr(prhs[2]);
    int npts = N[0];
    
    float u_max = 0;
	float v_max = 0;
	for (int i=0; i<npts; i++){
		if (p[2*i+0]>u_max) u_max=p[2*i+0];
		if (p[2*i+1]>v_max) v_max=p[2*i+1];
	}

	// allocate number of buckets needed
    int32_t bucket_width = dims[0];
    int32_t bucket_height= dims[1];
	int32_t bucket_cols = (int32_t)floor(u_max/bucket_width)  + 1;
	int32_t bucket_rows = (int32_t)floor(v_max/bucket_height) + 1;
    
	// assign matches to their buckets
    int *u = new int[npts];
    int *v = new int[npts];
	for (int i=0; i<npts; i++){
		u[i] = floor(p[2*i+0]/bucket_width);
		v[i] = floor(p[2*i+1]/bucket_height);
	}
    
    std::list<double> u_ind(u,u+npts);
    u_ind.sort();
    u_ind.unique();
    std::list<double> v_ind(u,u+npts);
    v_ind.sort();
    v_ind.unique();
    
    std::vector<std::vector<int>> b;
	b.resize(bucket_cols*bucket_rows);
    for (std::list<double>::iterator itu=u_ind.begin(); itu!=u_ind.end(); ++itu)
        for (std::list<double>::iterator itv=v_ind.begin(); itv!=v_ind.end(); ++itv)
            for (int i=0; i<npts; i++)
                if (u[i]==*itu & v[i]==*itv)
                    b[*itu * bucket_rows + *itv].push_back(i);
    
    //for (int i=0; i<b.size(); i++){
    //    for (int j=0; j<b[i].size(); j++)
    //        std::cout << b[i][j] << " ";
    //    std::cout << "\n";
    //}
    
    /* output cell to matlab */
    plhs[0] = mxCreateCellMatrix(b.size(), 1);
    mxArray *buckets = plhs[0];
    for(mwIndex i=0; i<b.size(); i++){
        mxArray *A;
        A = mxCreateDoubleMatrix(1, b[i].size(), mxREAL); // mxArray inside each Cell
        double *AData = mxGetPr(A);
        for (int j=0; j<b[i].size(); j++){
            AData[j] = b[i][j];
        }
        mxSetCell(buckets, i, A);
    }
    
    delete[] u, v;
}
