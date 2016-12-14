#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "mex.h"


/****** FUNCTIONS FOR DEBUGGING *********/

// Returns random double, uniform distribution in [-1,+1].
inline double frand() { return  (double)rand()/((int)1<<30)-1; }

// l = norml(l)
void mexFunction( int nargout,
                  mxArray **argout,
                  int nargin,
                  const mxArray **argin
                                   )
{
  if ( nargin!=1 ) mexErrMsgTxt("One input required.");
  const int M = mxGetM(argin[0]), N = mxGetN(argin[0]);
  if ( !M || !N ) {
    argout[0] = mxCreateDoubleMatrix(M,N,mxREAL);
    return;
  }
  if ( N!=3 ) mexErrMsgTxt("Input must be Mx3 matrix.");
  const double *l0 = (double*)mxGetPr(argin[0]), *l1 = l0+M, *l2 = l1+M;
  double *k0 = (double*)mxGetPr(argout[0] = mxCreateDoubleMatrix(M,3,mxREAL)), *k1 = k0+M, *k2 = k1+M;
  for ( int m=0; m<M; m++, l0++, l1++, l2++, k0++, k1++, k2++ ) {
    double r = sqrt(*l0**l0+*l1**l1);
    *k0 = *l0/r;
    *k1 = *l1/r;
    *k2 = *l2/r;
  }
}
