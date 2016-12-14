//Mexification of objective function for nonlinear minimization in line3d_from_lP_nonlin.m.

#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "mex.h"


/****** FUNCTIONS FOR DEBUGGING *********/

// Returns random double, uniform distribution in [-1,+1].
inline double frand() { return  (double)rand()/((int)1<<30)-1; }

void printm( const double x[][3], int N )
{
  printf("[");
  for ( int i=0; i<3; i++ ) {
    for ( int n=0; n<N; n++ ) printf("%e ",x[n][i]);
    printf("\n");
  }
  printf("]\n");
}

void printx( const double *x, int N )
{
  for ( int n=0; n<N; n++ ) printf("%.14g ",x[n]);
  printf("\n");
}

void printnormx( double x[3] ) 
{
  double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  printf("[%e %e %e]\n",x[0]/r,x[1]/r,x[2]/r);
}


/*********** MATRIX ALGEBRA ***********/

#define dot3(x,y)  ( (x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2] )
#define cross(x,y,r) \
  (r)[0] = (x)[1]*(y)[2] - (x)[2]*(y)[1]; \
  (r)[1] = (x)[2]*(y)[0] - (x)[0]*(y)[2]; \
  (r)[2] = (x)[0]*(y)[1] - (x)[1]*(y)[0];
#define  dot_op(s,x,y,r) \
  (r)[0] = (x)[0]s(y)[0]; \
  (r)[1] = (x)[1]s(y)[1]; \
  (r)[2] = (x)[2]s(y)[2];
#define reciproc(x,y) \
  (y)[0] = (x)[1]*(x)[2]; \
  (y)[1] = (x)[2]*(x)[0]; \
  (y)[2] = (x)[0]*(x)[1];
#define matrix_times(A,x,r) \
  (r)[0] = (A)[0][0]*(x)[0] + (A)[1][0]*(x)[1] + (A)[2][0]*(x)[2]; \
  (r)[1] = (A)[0][1]*(x)[0] + (A)[1][1]*(x)[1] + (A)[2][1]*(x)[2]; \
  (r)[2] = (A)[0][2]*(x)[0] + (A)[1][2]*(x)[1] + (A)[2][2]*(x)[2];
#define matrix_inv(A,B) \
  (B)[0][0] = (A)[1][1]*(A)[2][2] - (A)[1][2]*(A)[2][1]; \
  (B)[1][0] = (A)[1][2]*(A)[2][0] - (A)[1][0]*(A)[2][2]; \
  (B)[2][0] = (A)[1][0]*(A)[2][1] - (A)[1][1]*(A)[2][0]; \
  (B)[0][1] = (A)[2][1]*(A)[0][2] - (A)[2][2]*(A)[0][1]; \
  (B)[1][1] = (A)[2][2]*(A)[0][0] - (A)[2][0]*(A)[0][2]; \
  (B)[2][1] = (A)[2][0]*(A)[0][1] - (A)[2][1]*(A)[0][0]; \
  (B)[0][2] = (A)[0][1]*(A)[1][2] - (A)[0][2]*(A)[1][1]; \
  (B)[1][2] = (A)[0][2]*(A)[1][0] - (A)[0][0]*(A)[1][2]; \
  (B)[2][2] = (A)[0][0]*(A)[1][1] - (A)[0][1]*(A)[1][0];
#define xdif(x,y) \
  (y)[0] = (x)[2]-(x)[1]; \
  (y)[1] = (x)[0]-(x)[2]; \
  (y)[2] = (x)[1]-(x)[0]; \



void obj( double *y,
          const double *p,
          const double *Ppv,
          const double *AP,
          const double *u,
          const int K )
{
  // A1 = [1 p(1:2)']*AP(:,:,1);
  // A2 = [1 p(3:4)']*AP(:,:,2);
  // Lpv = [cross(A1(1:3)',A2(1:3)')' A2(4)*A1(1:3)-A1(4)*A2(1:3)];
  double
    A1[4] = { AP[ 0   ] + p[0]*AP[ 1   ] + p[1]*AP[ 2   ],
              AP[ 3   ] + p[0]*AP[ 4   ] + p[1]*AP[ 5   ],
              AP[ 6   ] + p[0]*AP[ 7   ] + p[1]*AP[ 8   ],
              AP[ 9   ] + p[0]*AP[10   ] + p[1]*AP[11   ] },
    A2[4] = { AP[ 0+12] + p[2]*AP[ 1+12] + p[3]*AP[ 2+12],
              AP[ 3+12] + p[2]*AP[ 4+12] + p[3]*AP[ 5+12],
              AP[ 6+12] + p[2]*AP[ 7+12] + p[3]*AP[ 8+12],
              AP[ 9+12] + p[2]*AP[10+12] + p[3]*AP[11+12] },
    Lpv[6];
  cross(A1,A2,Lpv);
  Lpv[3] = A2[3]*A1[0] - A1[3]*A2[0];
  Lpv[4] = A2[3]*A1[1] - A1[3]*A2[1];
  Lpv[5] = A2[3]*A1[2] - A1[3]*A2[2];

  //   y = zeros(3*K,1);
  //   for k = 1:K
  //     l = Lpv*Ppv(:,(1:3)+(k-1)*3);
  //     l = l/sqrt(l(1)^2+l(2)^2);
  //     y = [ y l*[u(:,k,1);1] l*[u(:,k,2);1] ];
  //   end
  for ( int k=0; k<K; k++ ) {
    const double *Pk = Ppv+18*k, *uk = u+2*k, *vk = uk+2*K;
    double l[3] = { Lpv[0]*Pk[0   ]+Lpv[1]*Pk[1   ]+Lpv[2]*Pk[2   ]+Lpv[3]*Pk[3   ]+Lpv[4]*Pk[4   ]+Lpv[5]*Pk[5   ],
                    Lpv[0]*Pk[0+ 6]+Lpv[1]*Pk[1+ 6]+Lpv[2]*Pk[2+ 6]+Lpv[3]*Pk[3+ 6]+Lpv[4]*Pk[4+ 6]+Lpv[5]*Pk[5+ 6],
                    Lpv[0]*Pk[0+12]+Lpv[1]*Pk[1+12]+Lpv[2]*Pk[2+12]+Lpv[3]*Pk[3+12]+Lpv[4]*Pk[4+12]+Lpv[5]*Pk[5+12] },
      rl = 1/sqrt(l[0]*l[0]+l[1]*l[1]),
      *yk = y+2*k;
    yk[0] = rl*( l[0]*uk[0] + l[1]*uk[1] + l[2] );
    yk[1] = rl*( l[0]*vk[0] + l[1]*vk[1] + l[2] );
  }
}


// [y,J] = obj_line3d_from_lP(p,Ppv,AP,cat(3,u,v));
void mexFunction( int nargout,
                  mxArray **argout,
                  int nargin,
                  const mxArray **argin
                                   )
{
  if ( nargin!=4 ) mexErrMsgTxt("4 inputs required.");
  const double
    *p = (double*)mxGetPr(argin[0]),
    *Ppv = (double*)mxGetPr(argin[1]),
    *AP = (double*)mxGetPr(argin[2]),
    *u = (double*)mxGetPr(argin[3]);
  const int *idims = mxGetDimensions(argin[3]), K = idims[1];  //printf("%i ",K);

  double *y = (double*)mxGetPr(argout[0] = mxCreateDoubleMatrix(2*K,1,mxREAL));
  obj(y,p,Ppv,AP,u,K);

  if ( nargout<2 ) return;

// dif = sqrt(eps);
// J = zeros(length(y),length(p));
// for i = 1:length(p)
//   pdif = p;
//   pdif(i) = pdif(i) + dif;
//   J(:,i) = (FFF(pdif,Ppv,AP,s) - y)/dif;
// end
  double *J = (double*)mxGetPr(argout[1] = mxCreateDoubleMatrix(2*K,4,mxREAL));
  for ( int i=0; i<4; i++ ) {
    double pdif[4] = { p[0],p[1],p[2],p[3] }, *Ji = J+2*K*i;
    pdif[i] += 1e-8;
    obj(Ji,pdif,Ppv,AP,u,K);
    for ( int j=0; j<2*K; j++ ) Ji[j] = (Ji[j]-y[j])*1e+8;
  }

}
