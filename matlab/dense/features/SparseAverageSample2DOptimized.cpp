/*
% *  This code was used in the following articles:
% *  [1] Learning 3-D Scene Structure from a Single Still Image, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on 3D Representation for Recognition (3dRR-07), 2007.
% *      (best paper)
% *  [2] 3-D Reconstruction from Sparse Views using Monocular Vision, 
% *      Ashutosh Saxena, Min Sun, Andrew Y. Ng, 
% *      In ICCV workshop on Virtual Representations and Modeling 
% *      of Large-scale environments (VRML), 2007. 
% *  [3] 3-D Depth Reconstruction from a Single Still Image, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      International Journal of Computer Vision (IJCV), Aug 2007. 
% *  [6] Learning Depth from Single Monocular Images, 
% *      Ashutosh Saxena, Sung H. Chung, Andrew Y. Ng. 
% *      In Neural Information Processing Systems (NIPS) 18, 2005.
% *
% *  These articles are available at:
% *  http://make3d.stanford.edu/publications
% * 
% *  We request that you cite the papers [1], [3] and [6] in any of
% *  your reports that uses this code. 
% *  Further, if you use the code in image3dstiching/ (multiple image version),
% *  then please cite [2].
% *  
% *  If you use the code in third_party/, then PLEASE CITE and follow the
% *  LICENSE OF THE CORRESPONDING THIRD PARTY CODE.
% *
% *  Finally, this code is for non-commercial use only.  For further 
% *  information and to obtain a copy of the license, see 
% *
% *  http://make3d.stanford.edu/publications/code
% *
% *  Also, the software distributed under the License is distributed on an 
% * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
% *  express or implied.   See the License for the specific language governing 
% *  permissions and limitations under the License.
% *
% */

#include <cstdlib>
#include <stdio.h>
#include <math.h>
/* #include <stdio.h> */
#include "mex.h"

/* Input Arguments */

#define	SMATRIX_IN	prhs[0]
#define	RowS	prhs[1]
#define	ColS	prhs[2]
#define	SAMPLE_P	prhs[3]
#define	SUP_MASK	prhs[4]


/* Output Arguments */

#define	SMATRIX_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif



/*    avgImg(SparseMOut,*rS,*cS,SamplePIn,SupMask,SparseMIn,m,n,Sm); */
static void avgImg(
		   double	SparseMOut[],
		   double	rS,
		   double	cS,
                   double       SamplePIn[],
                   double       SupMask[],
 		   double	SparseMIn[],
		   unsigned int m,
                   unsigned int n,
                   unsigned int Sm)
{
    double	sumOut;
    int CountNonZeros;
    int i,k,l;
    int i_adjusted, l_start, l_startIntom, l_end, innerCalc;
 
    /* Pre-calculate the fix number */
    int cSby2 = (int) floor(cS/2);
    int rSby2 = (int) floor(rS/2);
//    printf("cSby2 = %d", cSby2); //confirmed
 //   printf("rSby2 = %d", rSby2); // confirmed

    for (i=0;i<Sm;i++){ /*/ Count the No. of Sample points */
        sumOut=0;
        CountNonZeros = 0;
          
        /* Pre-calculate the fix number */
        int diffSamplePIn =(int) (SamplePIn[i+Sm] - 1) - cSby2; 
	l_start = MAX(diffSamplePIn, 0);
	l_startIntom = l_start*m;
	l_end = (int) MIN(diffSamplePIn + cS, n);
      
	diffSamplePIn = (int)(SamplePIn[i]-1) - rSby2;
 
            for (k= MAX(diffSamplePIn, 0);  k< (int)MIN(diffSamplePIn+rS,m);  k++){
                for (l=l_start,innerCalc=l_startIntom+k ;l<l_end; l++){

                    /* Pre-calculate the fix number */
    	            /*innerCalc = l*m+k;*/
                    if ( SupMask[innerCalc] ) {
                      sumOut += SparseMIn[innerCalc];
                      CountNonZeros++;
                    }
                    else{
//                      printf("SupMask 0\n");
                    }
    	            innerCalc += m;
                }
            }
          if (CountNonZeros == 0)
             {mexErrMsgTxt("CountNonZeros is zeros. It will cause error.");
             }  
            SparseMOut[i]=sumOut/CountNonZeros;/* /(rS*cS) */;
   }
    
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
    
{ 
    double *SparseMOut;
    double *SparseMIn; 
    double *SamplePIn; 
    double *SupMask; 
    double *rS,*cS; 
    unsigned int m,n,Supm,Supn,Sm,Sn; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 5) { 
	mexErrMsgTxt("Five input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /* Check the dimensions of the inputs. */ 
    
    m = mxGetM(RowS);
    n = mxGetN(RowS); 
    if (!mxIsDouble(RowS) || mxIsComplex(RowS) || 
	(MAX(m,n) != 1)) { 
	mexErrMsgTxt("The row size should be a scalar."); 
    } 
    m = mxGetM(ColS);
    n = mxGetN(ColS); 
    if (!mxIsDouble(ColS) || mxIsComplex(ColS) || 
	(MAX(m,n) != 1)) { 
	mexErrMsgTxt("The column size should be a scalar."); 
    } 
    m = mxGetM(SMATRIX_IN); 
    n = mxGetN(SMATRIX_IN);
    if (!mxIsDouble(SMATRIX_IN) || mxIsComplex(SMATRIX_IN)) { 
	mexErrMsgTxt("Input Sparse Matrix Format is not right."); 
    }
    Supm = mxGetM(SUP_MASK); 
    Supn = mxGetN(SUP_MASK);
    if (!mxIsDouble(SUP_MASK) || mxIsComplex(SUP_MASK) ||
        Supm !=m || Supn != n) { 
	mexErrMsgTxt("Input Superpixel Mask is not right."); 
    }
    Sm = mxGetM(SAMPLE_P); /*/Number of Rows equals to Number of points to calculate*/
    Sn = mxGetN(SAMPLE_P);
    if (!mxIsDouble(SAMPLE_P) || mxIsComplex(SAMPLE_P) ||
        Sn != 2) { 
	mexErrMsgTxt("Input Sample Point Index format should be Integer matrix of size M by 2.");
    }

    
    /* Create a matrix for the return argument */ 
    SMATRIX_OUT = mxCreateDoubleMatrix(Sm, 1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    SparseMOut = mxGetPr(SMATRIX_OUT);
    
    rS = mxGetPr(RowS); 
    cS = mxGetPr(ColS); 
    SparseMIn = mxGetPr(SMATRIX_IN);
    SamplePIn = mxGetPr(SAMPLE_P);
    SupMask = mxGetPr(SUP_MASK);
        
    /* Do the actual computations in a subroutine */
    avgImg(SparseMOut,*rS,*cS,SamplePIn,SupMask,SparseMIn,m,n,Sm); 
    return;
    
}


