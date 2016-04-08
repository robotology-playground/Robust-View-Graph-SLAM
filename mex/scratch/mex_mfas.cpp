
#include "mex.h"
#include "../src/mfas.cpp"
using namespace std;

/****************************/

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
    
    double *g;  /* 2xN */
    g = mxGetPr(prhs[0]);
    double *w;  /* 1xN */
    w = mxGetPr(prhs[1]);
    
    size_t ncols; /* matrix dimensions */
    ncols = mxGetN(prhs[0]);
    vector<Edge> edges (ncols, make_pair(-1, -1));
    vector<double> weight; weight.assign(ncols, 0);
    for (int i=0;i<ncols;i++){
        edges[i].first =g[2*i+0];
        edges[i].second=g[2*i+1];
        weight[i]=w[i];}
    
    /*************************************************************/
    map<int, int> reindexing_key;
    reindex_problem(edges, reindexing_key);
    
    flip_neg_edges(edges, weight);
                     
    vector<int> order;
    mfas_ratio(edges, weight, order);
    
    vector<double> broken;
    broken_weight(edges, weight, order, broken);
    
    /*************************************************************/
    plhs[0] = mxCreateDoubleMatrix(2, edges.size(), mxREAL);
    double *out0 = mxGetPr(plhs[0]);
    for (int i=0;i<weight.size();i++){
        out0[2*i+0] =edges[i].first;
        out0[2*i+1] =edges[i].second;}
    
    plhs[1] = mxCreateDoubleMatrix(1, weight.size(), mxREAL);
    double *out1 = mxGetPr(plhs[1]);
    for (int i=0;i<weight.size();i++){
        out1[i] =weight[i];}
    
    plhs[2] = mxCreateDoubleMatrix(1, order.size(), mxREAL);
    double *out2 = mxGetPr(plhs[2]);
    for (int i=0;i<order.size();i++){
        out2[i] =order[i];}
    
    plhs[3] = mxCreateDoubleMatrix(1, broken.size(), mxREAL);
    double *out3 = mxGetPr(plhs[3]);
    for (int i=0;i<broken.size();i++){
        out3[i] =broken[i];}
}