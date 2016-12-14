/* C = imcorrH(X,Y,H,W,u [,method]) */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

typedef int int2[2];
typedef int int3[3];
typedef double double2[2];
typedef double double3[3];


inline double BilinearInterpolation( double *X, int M, int N, double i, double j )
{
  int ii = (int)i, ji = (int)j, k = ii + ji*M;
  return
      (1-(j-=ji))*((1-(i-=ii))*X[k] +
                i *X[k+1]) +
         j *((1-i)*X[k+M] +
             i *X[k+1+M]);
}

inline double BilinearInterpolation_uchar( unsigned char *X, int M, int N, double i, double j )
{
  int ii = (int)i, ji = (int)j, k = ii + ji*M;
  return
      (1-(j-=ji))*((1-(i-=ii))*X[k] +
                i *X[k+1]) +
         j *((1-i)*X[k+M] +
             i *X[k+1+M]);
}


void Htimes( double *H, double2 u, double *x )
{
  x[0] = H[0]*u[0] + H[3]*u[1] + H[6];
  x[1] = H[1]*u[0] + H[4]*u[1] + H[7];
  x[2] = H[2]*u[0] + H[5]*u[1] + H[8];
}


// derivative A of homography v = nhom(H*hom(u)).
// Derivative of homography v = h(u) = nhom(H*[u;1]) is:
//   A = [eye(2) -v]/y(3)*H(:,1:2), where y = H*[u;1].
//
void dhom( double *H, // homography columwise
           double2 u, // point u
           double *A   // derivative, columnwise
  )
{
  double y[3];
  double2 v;

  Htimes(H,u,y);
  v[0] = y[0]/y[2];
  v[1] = y[1]/y[2];

  A[0] = ( H[0] - v[0]*H[2] )/y[2];
  A[1] = ( H[1] - v[1]*H[2] )/y[2];
  A[2] = ( H[3] - v[0]*H[5] )/y[2];
  A[3] = ( H[4] - v[0]*H[5] )/y[2];
}


double ncc( void *X, int2 NX, // image 1 and its size
            int Xuint8,       // flag for X being uint8
            void *Y, int2 NY, // image 2 and its size
            int Yuint8,       // flag for X being uint8
            double2 u,        // point in image 1 in which to compute ncc, in Matlab coord system
            double *H,         // homography matrix, columnwise, in Matlab coord system
            double *W, int NW, // displacement vectors/weights of window and their number
            int method         // 0 for nearest neighbor, 1 for bilin interpolation
   )
{
  #define INVALID -12345
  double A[4], EX, EY, DX, DY, XY, sumw, y[3];
  int i, nw;

  Htimes(H,u,y);
  double2 v = { y[0]/y[2], y[1]/y[2] }, // u1 transformed into Y
    d, e;

  // A := affinity which is derivative of homography in point u
  dhom(H,u,A);

  // compute ncc
  XY = EX = EY = DX = DY = sumw = 0;
  nw = 0;
  for ( i = 0; i < NW; i++ ) {
    double x, y,
           d1 = W[0+i*3], d2 = W[1+i*3], // displacement in image X
           w = W[2+3*i], // weight
           e1 = A[0]*d1 + A[2]*d2, e2 = A[1]*d1 + A[3]*d2, // displacement in image Y
           ud1 = u[0]+d1-1, ud2 = u[1]+d2-1,
           ve1 = v[0]+e1-1, ve2 = v[1]+e2-1; // final coordinates

    if ( ud1<1 || ud2<1 || ud1>=NX[0]-1 || ud2>=NX[1]-1 ||
         ve1<1 || ve2<1 || ve1>=NY[0]-1 || ve2>=NY[1]-1 ) continue;

    if ( method==1 ) {
      x = Xuint8 ?
        BilinearInterpolation_uchar((unsigned char*)X,NX[0],NX[1],ud1,ud2) :
        BilinearInterpolation      ((double       *)X,NX[0],NX[1],ud1,ud2);

      y = Yuint8 ?
        BilinearInterpolation_uchar((unsigned char*)Y,NY[0],NY[1],ve1,ve2) :
        BilinearInterpolation      ((double       *)Y,NY[0],NY[1],ve1,ve2);
    }
    else {
      x = Xuint8 ?
        ((unsigned char*)X)[(int)ud1+((int)ud2)*NX[0]] :
        ((double       *)X)[(int)ud1+((int)ud2)*NX[0]];

      y = Yuint8 ?
        ((unsigned char*)Y)[(int)ve1+((int)ve2)*NY[0]] :
        ((double       *)Y)[(int)ve1+((int)ve2)*NY[0]];
    }

    EX += x*w;   EY += y*w;
    DX += x*x*w; DY += y*y*w;
    XY += x*y*w;
    sumw += w;
    nw++;
    //printf("%f %f %f\n",x,y,w);
  }

  if ( nw > 2 ) { // for correlation we need at least 3 elements; for 2 elements it returned >1 or <-1
    EX /= sumw; EY /= sumw;
    DX = DX/sumw - EX*EX; DY = DY/sumw - EY*EY;
    XY = (XY/sumw - EX*EY) / sqrt(DX*DY);
    //printf("%f %f %i %f\n",DX,DY,nw,sumw);
    return ( (DX>1e-9) && (DY>1e-9) ) ? XY : mxGetNaN();
  }
  else
    return mxGetNaN();
}


void mexFunction( int nargout, mxArray *argout[], int nargin, const mxArray *argin[] )
{
  void *X, *Y;
  double *H, *W, *u, *C;
  int NW, Nu, i, Xuint8, Yuint8, method;
  int2 NX, NY;
  char method_str[20];

  #define  dst_C       ( argout[0] )
  #define  src_X      ( argin[0] )
  #define  src_Y      ( argin[1] )
  #define  src_H      ( argin[2] )
  #define  src_W      ( argin[3] )
  #define  src_u      ( argin[4] )
  #define  src_method  ( argin[5] )

  if ( (nargin!=5) && (nargin!=6) ) mexErrMsgTxt("Bad number of input parameters.");

  X = (void*)mxGetPr(src_X);
  NX[0] = mxGetM(src_X);
  NX[1] = mxGetN(src_X);
  switch ( mxGetClassID(src_X) ) {
    case mxDOUBLE_CLASS : Xuint8 = 0; break;
    case mxUINT8_CLASS : Xuint8 = 1; break;
    default : mexErrMsgTxt("X must be either double or uint8.");
  }

  Y = (void*)mxGetPr(src_Y);
  NY[0] = mxGetM(src_Y);
  NY[1] = mxGetN(src_Y);
  switch ( mxGetClassID(src_Y) ) {
    case mxDOUBLE_CLASS : Yuint8 = 0; break;
    case mxUINT8_CLASS : Yuint8 = 1; break;
    default : mexErrMsgTxt("Y must be either double or uint8.");
  }

  H = mxGetPr(src_H);
  if ( (mxGetM(src_H)!=3) || (mxGetN(src_H)!=3) ) mexErrMsgTxt("H must be 3x3 matrix.");

  W = mxGetPr(src_W);
  
  if ( mxGetM(src_W)!=3 ) mexErrMsgTxt("W must have 3 rows.");
  NW = mxGetN(src_W);
  
  u = mxGetPr(src_u);
  Nu = mxGetN(src_u); // Mu is checked below  

  if ( nargin==6 ) {
    mxGetString(src_method,method_str,19);
    method = strcmp(method_str,"nearest")==0 ? 0 : 1;
  }
  else method = 1;

  if ( (dst_C = mxCreateDoubleMatrix(1,Nu,mxREAL)) == 0 ) mexErrMsgTxt("Out of memory.");
  C = mxGetPr(dst_C);
  if ( mxIsEmpty(src_u) )
    return;
  if ( mxGetM(src_u)!=2 ) mexErrMsgTxt("Bad size of u.");

  for ( i = 0; i < Nu; i++ ) {
    double2 ui = { u[0+i*2], u[1+i*2] };
    C[i] = ncc(X,NX,Xuint8,Y,NY,Yuint8,ui,H,W,NW,method);
  }
  
}
