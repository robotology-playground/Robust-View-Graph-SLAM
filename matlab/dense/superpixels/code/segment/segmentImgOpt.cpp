/*Including the mexFunction as a matlab interface
usage:
a = segmentImgOpt( sigm, k, min, img, OltputFolder, PPMoption);
 Input--
 sigm
 k
 min
 img: rgb 3 color channel 2d matrix
 OltputFolder : path of the  output folder
 PPMoption : 1 if storaging PPM file
 
 Output--
 a : unsorted index superpixel 2d matrix*/



#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "segment-image-opt.h"
#include "mex.h"
//#include <cstdio>

/* */
#define PPM_PATH        prhs[4]
#define PPM_OPTION        prhs[5]

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  if(nrhs!=6) {
    mexErrMsgTxt("usage: segmentedImage = SEGMENT(sigma, k, min, inputimage, OltputFolder, PPMoption)");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments");
  }
  if (!mxIsChar( PPM_PATH) || mxIsComplex( PPM_PATH) ) {
        mexErrMsgTxt("WRL path is not a string.");
    }

  const int *dims;  
  int height, width, ndim;
  char *PpmFolder;
  float sigma = (float)*mxGetPr(prhs[0]);
  float k = (float)*mxGetPr(prhs[1]);
  int min_size = (int)*mxGetPr(prhs[2]);
  unsigned char *Imptr;
  double *SegOut;
  PpmFolder = mxArrayToString(PPM_PATH);
  double *PpmOption = mxGetPr( PPM_OPTION);
/*  float *Imptr;*/
  unsigned short int *output;
  
  height=(mxGetDimensions(prhs[3]))[0];
  width=(mxGetDimensions(prhs[3]))[1];

  Imptr = (unsigned char*) mxGetPr(prhs[3]);
/*  Imptr = (float*) mxGetPr(prhs[3]);*/
/*  image<rgb> *im = new image<rgb>(width, height, true);*/
/*  image<unsigned char> *r = new image<unsigned char>(width, height, true);
  image<unsigned char> *g = new image<unsigned char>(width, height, true);
  image<unsigned char> *b = new image<unsigned char>(width, height, true);*/

/*  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {*/
/*      imRef(r, x, y) = Imptr[height*x+y];
      imRef(g, x, y) = Imptr[height*x+y];
      imRef(b, x, y) = Imptr[height*x+y];*/
/*      imRef(g, x, y) = Imptr[height*x+y + (width* height * 1 -1)];
      imRef(b, x, y) = Imptr[height*x+y + (width* height * 2 -1)];*/
/*      imRef(r, x, y) = Imptr[height*x+y];
      imRef(g, x, y) = Imptr[height*x+y + (width* height * 1 -1)];
      imRef(b, x, y) = Imptr[height*x+y + (width* height * 2 -1)];*/
/*      imRef(im, x, y).r = Imptr[height*x+y];
      imRef(im, x, y).g = Imptr[height*x+y + (width* height * 1 -1)];
      imRef(im, x, y).b = Imptr[height*x+y + (width* height * 2 -1)];*/
/*    }
  }*/
/*  mexWarnMsgTxt(" r g b setup OK");*/
/*  memcpy(im->data, mxGetData(prhs[3]), height* width * sizeof(rgb));*/
  plhs[0] = mxCreateDoubleMatrix ( height, width, mxREAL);
  SegOut =  mxGetPr(plhs[0]);

/*  image<rgb> *seg = segment_image( im, SegOut, height, width, sigma, k, min_size);*/
  image<rgb> *seg = segment_image( Imptr, SegOut, height, width, sigma, k, min_size, PpmOption);

  if (*PpmOption){
      savePPM( seg, PpmFolder);
  }

  /*delete [] r;
  delete [] g;
  delete [] b;*/
  return;
}
