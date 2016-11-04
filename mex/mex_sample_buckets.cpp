#include "mex.h"
#include <iostream> // std::cout
#include <string>
#include <math.h>       /* floor */
#include <list>
#include <vector>
#include <memory.h>
#include <algorithm>    // std::random_shuffle

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    /* Convert inputs into standard C++ format without using mex.h */
    mxArray *bucket;
    double *ptr;
    const mwSize *dims = mxGetDimensions(prhs[0]);
    std::vector<std::vector<int>> b;
    b.resize(dims[0]);
    for (int jcell=0; jcell<dims[0]; jcell++){
        //mexPrintf("Convering backet %d \n", jcell);
        bucket = mxGetCell(prhs[0],jcell);
        const mwSize *dims2 = mxGetDimensions(bucket);
        int number_of_dims = mxGetNumberOfDimensions(bucket);
        ptr = mxGetPr(bucket);
        for (int c=0; c<dims2[1]; c++){
            //mexPrintf("Getting index %d \n", c);
            b[jcell].push_back(*ptr++);
        }
    }
    
    /* Generate samples from buckets */
    int buck = 0;
    int tri = 0;
    std::vector<int> ind;
    while (ind.size()<std::max(5, (int)b.size()) & tri<20){
        tri++;
        if (buck>=b.size()){ // checks if end of buckets is reached.
            if (ind.size()>4) // finished all buckets and have at least 5 points? Then, exit.
                break;
            buck=1+buck%2; // finished with less than 5 points? Then, repeat bucketing using even buckets
        }
        //mexPrintf("Processing backet %d with %d points \n", buck, b[buck].size());
        if (b[buck].size()<3){ // bucket is empty? Then, move to next one
            buck+=2;
            continue;
        }
        int idx = 0;
        std::random_shuffle(b[buck].begin(),b[buck].end());
        while (idx<b[buck].size()){
            int k=0;
            for (int i=0; i<ind.size(); i++){ // check if point index is new
                if (ind[i]==b[buck][idx])
                    k++;
            }
            if (k<1)
                ind.push_back(b[buck][idx]);
            idx++;
        }
        buck+=2;
    }
    
    /* Convert output into mex format */
    plhs[0] = mxCreateDoubleMatrix(1,ind.size(),mxREAL);
    double *samples = mxGetPr(plhs[0]); /* output samples */
    std::random_shuffle(ind.begin(),ind.end());
    for (int i=0; i<ind.size(); i++){
        samples[i] = ind[i];
    }

}