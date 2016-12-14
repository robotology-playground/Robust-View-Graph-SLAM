#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "mex.h"


/*********** MATRIX ALGEBRA ***********/

#define cross(x,y,r) \
  (r)[0] = (x)[1]*(y)[2] - (x)[2]*(y)[1]; \
  (r)[1] = (x)[2]*(y)[0] - (x)[0]*(y)[2]; \
  (r)[2] = (x)[0]*(y)[1] - (x)[1]*(y)[0];


// c = cros(a)
// c = cros(a,b)
void mexFunction( int nargout,
                  mxArray **argout,
                  int nargin,
                  const mxArray **argin
                                   )
{
  if ( nargin==1 ) {
    const int M = mxGetM(argin[0]), N = mxGetN(argin[0]);
    if ( !M || !N ) {
      argout[0] = mxCreateDoubleMatrix(M,N,mxREAL);
      return;
    }
    if ( M!=3 || N!=2 ) mexErrMsgTxt("Input must be 3x2 vector.");
    const double *a = (double*)mxGetPr(argin[0]);
    double *c = (double*)mxGetPr(argout[0] = mxCreateDoubleMatrix(3,1,mxREAL));
    cross(a,a+3,c);
  }
  else if ( nargin==2 ) {
    const int M = mxGetM(argin[0]), N = mxGetN(argin[0]);
    if ( !M || !N ) {
      argout[0] = mxCreateDoubleMatrix(M,N,mxREAL);
      return;
    }
    if ( M!=3 || mxGetM(argin[1])!=3 || mxGetN(argin[1])!=N ) mexErrMsgTxt("Wrong size of input vectors.");
    const double *a = (double*)mxGetPr(argin[0]),
                 *b = (double*)mxGetPr(argin[1]);
    double *c = (double*)mxGetPr(argout[0] = mxCreateDoubleMatrix(3,N,mxREAL));
    for ( int n=0; n<N; n++, a+=3, b+=3, c+=3 ) { cross(a,b,c); }
  }
  else mexErrMsgTxt("One or two input arguments required.");
}
