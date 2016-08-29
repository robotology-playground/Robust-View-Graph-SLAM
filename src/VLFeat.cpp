/** @internal
** @file     sift.c
** @author   Andrea Vedaldi
** @brief    Scale Invariant Feature Transform (SIFT) - MEX
**/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

//extern "C" {
//#include <mexutils.h>
#include <vl/mathop.h>
#include "vl/imopv.h"
#include <vl/sift.h>

#include <math.h>
#include <assert.h>

#include <cxcore.h>
#include <highgui.h>

#include <iostream>

//#include "opencv2/calib3d.hpp"

#include "epnp.h"

#include <vl/generic.h>

#include<stdlib.h>
#include<string.h>
//}

typedef struct
{
	int k1 ;
	int k2 ;
	double score ;
} Pair ;


VL_INLINE void
	transpose_descriptor (vl_sift_pix* dst, vl_sift_pix* src)
{
	int const BO = 8 ;  /* number of orientation bins */
	int const BP = 4 ;  /* number of spatial bins     */
	int i, j, t ;

	for (j = 0 ; j < BP ; ++j) {
		int jp = BP - 1 - j ;
		for (i = 0 ; i < BP ; ++i) {
			int o  = BO * i + BP*BO * j  ;
			int op = BO * i + BP*BO * jp ;
			dst [op] = src[o] ;
			for (t = 1 ; t < BO ; ++t)
				dst [BO - t + op] = src [t + o] ;
		}
	}
}


static int
	korder (void const* a, void const* b) {
		double x = ((double*) a) [2] - ((double*) b) [2] ;
		if (x < 0) return -1 ;
		if (x > 0) return +1 ;
		return 0 ;
}


vl_bool
	check_sorted (double const * keys, vl_size nkeys)
{
	vl_uindex k ;
	for (k = 0 ; k + 1 < nkeys ; ++ k) {
		if (korder(keys, keys + 4) > 0) {
			return VL_FALSE ;
		}
		keys += 4 ;
	}
	return VL_TRUE ;
}


void VLSIFT(cv::Mat image, uint8_t* DATAdescr, double* DATAframes, int* nframes){
	//Take IplImage -> convert to SINGLE (float):
	float* frame = (float*)malloc(image.rows*image.cols*sizeof(float));
	uchar* Ldata      = (uchar *)image.data;
	for(int i=0;i<image.rows;i++)
		for(int j=0;j<image.cols;j++)
			frame[j*image.rows+i] = (float)Ldata[j*image.rows+i];
			
	/*
	FILE *fpp = fopen("c:\\Picture.txt", "w");
		for(int p=0;p<image->height*image->width; p++){
			fprintf(fpp, "%f\n",frame[p] );
		}
		fclose(fpp);
		*/

	// VL SIFT computation:
	vl_sift_pix const *data ;
	int M, N ;
	data = (vl_sift_pix*)frame;
	M = image.rows;
	N = image.cols;

	int verbose = 1 ;
	int O     =   -1 ; //Octaves
	int S     =   3 ; //Levels
	int                o_min =   0 ;

	double             edge_thresh = -1 ;
	double             peak_thresh =  -1 ;
	double             norm_thresh = -1 ;
	double             magnif      = -1 ;
	double             window_size = -1 ;

	//mxArray           *ikeys_array = 0 ; //?
	double            *ikeys = 0 ; //?
	int                nikeys = -1 ; //?
	vl_bool            force_orientations = 0 ;
	vl_bool            floatDescriptors = 0 ;

	/* -----------------------------------------------------------------
	*                                                            Do job
	* -------------------------------------------------------------- */
	{
		VlSiftFilt        *filt ;
		vl_bool            first ;
		double            *frames = 0 ;
		uint8_t              *descr  = 0 ;
		int                reserved = 0, i,j,q ;

      //std::cout << "1" << std::endl;

		/* create a filter to process the image */
		filt = vl_sift_new (M, N, O, S, o_min) ;
		
		//std::cout << "2" << std::endl;

		if (peak_thresh >= 0) vl_sift_set_peak_thresh (filt, peak_thresh) ;
		if (edge_thresh >= 0) vl_sift_set_edge_thresh (filt, edge_thresh) ;
		if (norm_thresh >= 0) vl_sift_set_norm_thresh (filt, norm_thresh) ;
		if (magnif      >= 0) vl_sift_set_magnif      (filt, magnif) ;
		if (window_size >= 0) vl_sift_set_window_size (filt, window_size) ;
		
			

		//if (verbose) {
		//	printf("vl_sift: filter settings:\n") ;
		//	printf("vl_sift:   octaves      (O)      = %d\n",
		//		vl_sift_get_noctaves      (filt)) ;
		//	printf("vl_sift:   levels       (S)      = %d\n",
		//		vl_sift_get_nlevels       (filt)) ;
		//	printf("vl_sift:   first octave (o_min)  = %d\n",
		//		vl_sift_get_octave_first  (filt)) ;
		//	printf("vl_sift:   edge thresh           = %g\n",
		//		vl_sift_get_edge_thresh   (filt)) ;
		//	printf("vl_sift:   peak thresh           = %g\n",
		//		vl_sift_get_peak_thresh   (filt)) ;
		//	printf("vl_sift:   norm thresh           = %g\n",
		//		vl_sift_get_norm_thresh   (filt)) ;
		//	printf("vl_sift:   window size           = %g\n",
		//		vl_sift_get_window_size   (filt)) ;
		//	printf("vl_sift:   float descriptor      = %d\n",
		//		floatDescriptors) ;

		//	printf((nikeys >= 0) ?
		//		"vl_sift: will source frames? yes (%d read)\n" :
		//	"vl_sift: will source frames? no\n", nikeys) ;
		//	printf("vl_sift: will force orientations? %s\n",
		//		force_orientations ? "yes" : "no") ;
		//}

		/* ...............................................................
		*                                             Process each octave
		* ............................................................ */
		i     = 0 ;
		first = 1 ;
		while (1) {
			int                   err ;
			VlSiftKeypoint const *keys  = 0 ;
			int                   nkeys = 0 ;

			//if (verbose) {
			//	printf ("vl_sift: processing octave %d\n",
			//		vl_sift_get_octave_index (filt)) ;
			//}

			/* Calculate the GSS for the next octave .................... */
			if (first) {
				err   = vl_sift_process_first_octave (filt, data) ;
				first = 0 ;
			} else {
				err   = vl_sift_process_next_octave  (filt) ;
			}

			if (err) break ;

		//	if (verbose > 1) {
		//		printf("vl_sift: GSS octave %d computed\n",
		//			vl_sift_get_octave_index (filt));
		//	}

			/* Run detector ............................................. */
			if (nikeys < 0) {
				vl_sift_detect (filt) ;

				keys  = vl_sift_get_keypoints  (filt) ;
				nkeys = vl_sift_get_nkeypoints (filt) ;
				i     = 0 ;

			//	if (verbose > 1) {
			//		printf ("vl_sift: detected %d (unoriented) keypoints\n", nkeys) ;
			//	}
			} else {
				nkeys = nikeys ;
			}

			/* For each keypoint ........................................ */
			for (; i < nkeys ; ++i) {
				double                angles [4] ;
				int                   nangles ;
				VlSiftKeypoint        ik ;
				VlSiftKeypoint const *k ;

				/* Obtain keypoint orientations ........................... */
				if (nikeys >= 0) {
					vl_sift_keypoint_init (filt, &ik,
						ikeys [4 * i + 1] - 1,
						ikeys [4 * i + 0] - 1,
						ikeys [4 * i + 2]) ;

					if (ik.o != vl_sift_get_octave_index (filt)) {
						break ;
					}

					k = &ik ;

					/* optionally compute orientations too */
					if (force_orientations) {
						nangles = vl_sift_calc_keypoint_orientations
							(filt, angles, k) ;
					} else {
						angles [0] = VL_PI / 2 - ikeys [4 * i + 3] ;
						nangles    = 1 ;
					}
				} else {
					k = keys + i ;
					nangles = vl_sift_calc_keypoint_orientations
						(filt, angles, k) ;
				}

				/* For each orientation ................................... */
				for (q = 0 ; q < nangles ; ++q) {
					vl_sift_pix  buf [128] ;
					vl_sift_pix rbuf [128] ;

					/* compute descriptor (if necessary) */
					vl_sift_calc_keypoint_descriptor (filt, buf, k, angles [q]) ;
					transpose_descriptor (rbuf, buf) ;

					/* make enough room for all these keypoints and more */
					if (reserved < (*nframes) + 1) {
						reserved += 2 * nkeys ;
						frames = (double*)realloc (frames, 4 * sizeof(double) * reserved) ;
						descr  = (uint8_t*)realloc (descr,  128 * sizeof(uint8_t) * reserved) ;
					}

					/* Save back with MATLAB conventions. Notice tha the input
					* image was the transpose of the actual image. */
					frames [4 * (*nframes) + 0] = k -> y ;
					frames [4 * (*nframes) + 1] = k -> x ;
					frames [4 * (*nframes) + 2] = k -> sigma ;
					frames [4 * (*nframes) + 3] = VL_PI / 2 - angles [q] ;


					for (j = 0 ; j < 128 ; ++j) {
						float x = 512.0F * rbuf [j] ;
						x = (x < 255.0F) ? x : 255.0F ;
						descr[128 * (*nframes) + j] = (uint8_t)x ;
					}

					++ (*nframes) ;
				} /* next orientation */
			} /* next keypoint */
		} /* next octave */

		//std::cout << "3" << std::endl;

	//	if (verbose) {
	//		printf ("vl_sift: found %d keypoints\n", (*nframes)) ;
	//	}
		// save variables:
		memcpy(DATAframes, frames, 4 * (*nframes ) * sizeof(double));
		memcpy(DATAdescr, descr, 128 * (*nframes ) * sizeof(uint8_t));
		
		/*
		FILE *fpd = fopen("c:\\Descr.txt", "w");
		for(int p=0;p<(*nframes)*128; p++){
			fprintf(fpd, "%f\n",(double)descr[p] );
		}
		fclose(fpd);

		FILE *fpf = fopen("c:\\Frames.txt", "w");
		for(int p=0;p<(*nframes)*4; p++){
			fprintf(fpf, "%f\n",frames[p] );
		}
		fclose(fpf);
		*/

		/* cleanup */
		vl_sift_delete (filt) ;
	} /* end: do job */


	return;
}

void VLMATCH(uint8_t* L1_pt,uint8_t* L2_pt, int K1, int K2, double thresh, int* nMatches, double* MATCHES ){
	//Match descriptors!

	//double thresh ;
	//int  K1, K2, ND ;
	//K1 = Tnframes;
	//K2 = Qnframes;
	int ND = 128;
	//uint8_t* L1_pt  ;
	//uint8_t* L2_pt ;
	//int* nMatches;
	//double* MATCHES;

	Pair* pairs_begin = (Pair*) malloc(sizeof(Pair) * (K1+K2)) ;
	Pair* pairs_iterator = pairs_begin ;

	int k1, k2 ;                                                        
	const int maxval = 0x7fffffff ;                         
	for(k1 = 0 ; k1 < K1 ; ++k1, L1_pt += ND ) {    //kalooo!                    

		int best = maxval ;                                     
		int second_best = maxval ;                              
		int bestk = -1 ;                                                  

		/* For each point P2[k2] in the second image... */                
		for(k2 =  0 ; k2 < K2 ; ++k2, L2_pt += ND) {                      

			int bin ;                                                       
			int acc = 0 ;                                         
			for(bin = 0 ; bin < ND ; ++bin) {                               
				int delta =                                         
					((int) L1_pt[bin]) -                              
					((int) L2_pt[bin]) ;                              
				acc += delta*delta ;                                          
			}                                                               

			/* Filter the best and second best matching point. */           
			if(acc < best) {                                                
				second_best = best ;                                          
				best = acc ;                                                  
				bestk = k2 ;                                                  
			} else if(acc < second_best) {                                  
				second_best = acc ;                                           
			}                                                               
		}                                                                 

		L2_pt -= ND*K2 ;                                                  

		/* Lowe's method: accept the match only if unique. */             
		if(thresh * (float) best < (float) second_best &&                 
			bestk != -1) {                                                 
				pairs_iterator->k1 = k1 ;                                       
				pairs_iterator->k2 = bestk ;                                    
				pairs_iterator->score = best ;                                  
				pairs_iterator++ ;  
				(*nMatches)++;
		}                                                                 
	}                                                                   

	Pair* pairs_end = pairs_iterator ;
	//double* M_pt = (double*)calloc((pairs_end-pairs_begin)*2,sizeof(double));
	double* M_pt = (double*)calloc((*nMatches)*2,sizeof(double));
	//double* M_start = M_pt;

	for(pairs_iterator = pairs_begin ;
		pairs_iterator < pairs_end  ;
		++pairs_iterator) {
			*M_pt++ = pairs_iterator->k1 ;
			*M_pt++ = pairs_iterator->k2 ;
	}
	M_pt -= (*nMatches)*2 ;
	memcpy(MATCHES,M_pt,(*nMatches) * 2 * sizeof(double));
	free(pairs_begin) ;
	free(M_pt);

	return;
}

typedef unsigned char uchar;
//int main()
//{
	// Load template image:
//	IplImage* Timage = cvLoadImage("T.png",0);

//	double*            TFrames = (double*)calloc ( 4 * 10000, sizeof(double) ) ;
//	uint8_t*          TDescr  = (uint8_t*)calloc ( 128 * 10000, sizeof(uint8_t) ) ;
//	int                Tnframes = 0;
//	VLSIFT(Timage, TDescr, TFrames, &Tnframes);
//	TFrames = (double*)realloc (TFrames, 4 * sizeof(double) * Tnframes) ; // = Y X Scale Angle
//	TDescr  = (uint8_t*)realloc (TDescr,  128 * sizeof(uint8_t) * Tnframes) ;
	
	/*
	for(int i=0;i<Tnframes;i++){
		cvCircle(Timage,                
		cvPoint(TFrames[0+i*4], TFrames[1+i*4]), TFrames[2+i*4],   
		cvScalar(255, 0, 0, 0),  
		1, 8, 0);  
	}
	cvShowImage("FrameT", Timage);
	cvWaitKey(0);
	*/


	// Load query frame:
//	IplImage* Qimage = cvLoadImage("Q2.png",0);

//	double*            QFrames = (double*)calloc ( 4 * 10000, sizeof(double) ) ;
//	uint8_t*          QDescr  = (uint8_t*)calloc ( 128 * 10000, sizeof(uint8_t) ) ;
//	int                Qnframes = 0;
//	VLSIFT(Qimage, QDescr, QFrames, &Qnframes);
//	QFrames = (double*)realloc (QFrames, 4 * sizeof(double) * Qnframes) ;
//	QDescr  = (uint8_t*)realloc (QDescr,  128 * sizeof(uint8_t) * Qnframes) ;
	
	/*
	for(int i=0;i<Qnframes;i++){
	cvCircle(Qimage,                
	cvPoint(QFrames[0+i*4], QFrames[1+i*4]), QFrames[2+i*4],   
	cvScalar(255, 0, 0, 0),  
	1, 8, 0);  
	}
	cvShowImage("FrameQ", Qimage);
	cvWaitKey(0);
	cvDestroyWindow("FrameQ");
	cvDestroyWindow("FrameT");
	*/

	// Match process:
//	int matchesFound=0;
//	double* MATCHES = (double*)calloc( 10000, sizeof(double) ) ;
//	VLMATCH(TDescr,QDescr, Tnframes, Qnframes, 5, &matchesFound, MATCHES );
//	MATCHES = (double*)realloc(MATCHES, sizeof(double) * matchesFound * 2) ;

	//Display matches for user checking:
	/*
	for(int i=0;i<Tnframes;i++){
		for(int m=0;m<matchesFound;m++){
			if(i == MATCHES[m*2]){
				cvCircle(Timage,                
					cvPoint(TFrames[0+i*4], TFrames[1+i*4]), TFrames[2+i*4],   
					cvScalar(255, 0, 0, 0),  
					1, 8, 0); 
			}
		}
	}
	cvShowImage("FrameT", Timage);

	for(int i=0;i<Qnframes;i++){
		for(int m=0;m<matchesFound;m++){
			if(i == MATCHES[m*2+1]){
				cvCircle(Qimage,                
					cvPoint(QFrames[0+i*4], QFrames[1+i*4]), QFrames[2+i*4],   
					cvScalar(255, 0, 0, 0),  
					1, 8, 0); 
			}
		}
	}
	cvShowImage("FrameQ", Qimage);
	cvWaitKey(0);
	cvDestroyWindow("FrameQ");
	cvDestroyWindow("FrameT");
	*/
//	getchar();
//	return(0);
//}
