// g = vgg_conv_lineseg(I,h,u1,u2)
//   (c) T. Werner, Nov 2001

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>


typedef int int2[2];
typedef double double2[2];
#define fround(x) ( (x)>0 ? (x)+.5 : (x)-.5 )


inline double BilinearInterpolation( double *X, int M, int N, double i, double j, double xInvalid )
{
  int ii = (int)i, ji = (int)j, k = ii + ji*M;

  return (ii>=0 && ji>=0 && ii<M-1 && ji<N-1) ?
      (1-(j-=ji))*((1-(i-=ii))*X[k] +
                i *X[k+1]) +
         j *((1-i)*X[k+M] +
                i *X[k+1+M]) :
      xInvalid;
}


inline double BilinearInterpolation_uchar( unsigned char *X, int M, int N, double i, double j, double xInvalid )
{
  int ii = (int)i, ji = (int)j, k = ii + ji*M;

  return (ii>=0 && ji>=0 && ii<M-1 && ji<N-1) ?
      (1-(j-=ji))*((1-(i-=ii))*X[k] +
                i *X[k+1]) +
         j *((1-i)*X[k+M] +
                i *X[k+1+M]) :
      xInvalid;
}


/* Image gradient along a line segment */
double imgraduv(
   void *I, int2 NI, // image and its size
   double *h, int rh, // derivative of gausian kernel, its size is 2*rh+1
   double2 v1, double2 v2, // line segment end points
   int UCHAR // says whether I is double or unsigned char
)
{
   double g, sumi, len, f, hi;
   double2 s, v;
   int nsumi, i, j, lenint;

   s[0] = v2[0] - v1[0];
   s[1] = v2[1] - v1[1];
   len = sqrt(s[0]*s[0] + s[1]*s[1]);
   lenint = int(fround(len));
   s[0] /= len;
   s[1] /= len;
//      printf("%g %g, %g %g, %g %g, %g\n",v1[0],v1[1],v2[0],v2[1],s[0],s[1],len);
   g = 0;
   for ( i = -rh; i <= rh; i++ ) {
     hi = h[i+rh];
     sumi = 0;
     nsumi = 0;
     for ( j = 0; j < lenint; j++ ) {
       v[0] = -s[1]*i + s[0]*j + v1[0];
       v[1] =  s[0]*i + s[1]*j + v1[1];
//	    printf("%i %i, %g %g\n",i,j,v[0],v[1]);
       f = UCHAR ?
         BilinearInterpolation_uchar((unsigned char*)I,NI[0],NI[1],v[0],v[1],-1234567) :
         BilinearInterpolation((double*)I,NI[0],NI[1],v[0],v[1],-1234567);
       if ( f != -1234567 ) {
          sumi += hi*f;
	  nsumi++;
       }
     }
     g += sumi/nsumi;
   }
   return g;
}


void mexFunction( int nargout, mxArray *argout[], int nargin, const mxArray *argin[] )
{
   void *I;
   double *h, *u1, *u2, *g;
   int rh, Nu, n, UCHAR;
   int2 NI;
   double2 v1, v2;

   if ( nargin != 4 ) mexErrMsgTxt("Bad number of input arguments.");

   if ( mxGetNumberOfDimensions(argin[0])!=2 ) mexErrMsgTxt("I must be 2D array.");
   NI[0] = mxGetM(argin[0]);
   NI[1] = mxGetN(argin[0]);

   switch ( mxGetClassID(argin[0]) ) {
   case mxDOUBLE_CLASS : UCHAR = 0; break;
   case mxUINT8_CLASS : UCHAR = 1; break;
   default : mexErrMsgTxt("I must be either double or uint8.");
   }

   I = mxGetPr(argin[0]);
   if ( !mxIsDouble(argin[1]) ) mexErrMsgTxt("h must be double vector.");
   rh = mxGetNumberOfElements(argin[1]);
   if ( rh%2==0 ) mexErrMsgTxt("Vector h must have odd length.");
   rh = (rh-1)/2;
   h = mxGetPr(argin[1]);
   if ( mxGetNumberOfDimensions(argin[2])!=2 || !mxIsDouble(argin[2]) || mxGetM(argin[2])!=2 ||
        mxGetNumberOfDimensions(argin[3])!=2 || !mxIsDouble(argin[3]) || mxGetM(argin[3])!=2 ||
        mxGetN(argin[2])!=mxGetN(argin[3]) ) mexErrMsgTxt("Bad size(s) of u1, u2.");
   Nu = mxGetN(argin[2]);
   u1 = mxGetPr(argin[2]);
   u2 = mxGetPr(argin[3]);

   argout[0] = mxCreateDoubleMatrix(1,Nu,mxREAL);
   if ( argout[0] == 0 ) mexErrMsgTxt("Out of memory.");
   g = mxGetPr(argout[0]);

   for ( n = 0; n < Nu; n++ ) {
      double2 v1, v2;
      v1[0] = u1[0+2*n] - 1;
      v1[1] = u1[1+2*n] - 1;
      v2[0] = u2[0+2*n] - 1;
      v2[1] = u2[1+2*n] - 1;
      g[n] = imgraduv(I,NI,h,rh,v1,v2,UCHAR);
   }

}
