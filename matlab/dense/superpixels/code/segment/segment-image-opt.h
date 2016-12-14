/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "filter.h"
#include "segment-graph.h"

// random color
rgb random_rgb(){ 
  rgb c;
  double r;
  
  c.r = (uchar)random();
  c.g = (uchar)random();
  c.b = (uchar)random();

  return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
			 int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
	      square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
	      square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
/*image<rgb> *segment_image(image<rgb> *im, int SegOut[], int height, int width, float sigma, float c, int min_size,
			  int *num_ccs) {*/
/*void segment_image(image<rgb> *im, double SegOut[], int height, int width, float sigma, float c, int min_size
			  ) {*/
/*void segment_image( image<unsigned char> *r, image<unsigned char> *g, image<unsigned char> *b, double SegOut[], int height, int width, float sigma, float c, int min_size
			 ) {*/
image<rgb> * segment_image( unsigned char *Imptr, double SegOut[], int height, int width, float sigma, float c, int min_size, double *PpmOption
			 ) {
/*image<rgb> *segment_image( image<float> *r, image<float> *g, image<float> *b, int height, int width, float sigma, float c, int min_size,*/
/*  int width = im->width();
  int height = im->height();*/

  image<float> *r = new image<float>(width, height);
  image<float> *g = new image<float>(width, height);
  image<float> *b = new image<float>(width, height);

  // smooth each color channel  
  int Area = width* height;
  for (int x = 0; x < width; x++) {
    int HightTimesX = height*x;
    for (int y = 0; y < height; y++) {
      int temp = HightTimesX+y;
      imRef(r, x, y) = Imptr[temp];
      temp += Area;
      imRef(g, x, y) = Imptr[temp];
      temp += Area;
      imRef(b, x, y) = Imptr[temp];

/*      imRef(g, x, y) = Imptr[temp + (Area * 1 -1)];;
      imRef(b, x, y) = Imptr[temp + (Area * 2 -1)];;*/
/*      imRef(r, x, y) = imRef(im, x, y).r;
      imRef(g, x, y) = imRef(im, x, y).g;
      imRef(b, x, y) = imRef(im, x, y).b;*/
    }
  }
  image<float> *smooth_r = smooth(r, sigma);
  image<float> *smooth_g = smooth(g, sigma);
  image<float> *smooth_b = smooth(b, sigma);
  delete r;
  delete g;
  delete b;
 
  // build graph
  edge *edges = new edge[width*height*4];
  int num = 0;
  for (int y = 0; y < height; y++) {
    int YtimesWidth = y * width;
    for (int x = 0; x < width; x++) {
      if (x < width-1) {
	edges[num].a = YtimesWidth + x;
	edges[num].b = YtimesWidth + (x+1);
	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	num++;
      }

      if (y < height-1) {
	edges[num].a = YtimesWidth + x;
	edges[num].b = YtimesWidth + width + x;
	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
	num++;
      }

      if ((x < width-1) && (y < height-1)) {
	edges[num].a = YtimesWidth + x;
	edges[num].b = YtimesWidth + width + (x+1);
	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
	num++;
      }

      if ((x < width-1) && (y > 0)) {
	edges[num].a = YtimesWidth + x;
	edges[num].b = YtimesWidth - width + (x+1);
	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
	num++;
      }
    }
  }
  delete smooth_r;
  delete smooth_g;
  delete smooth_b;

  // segment
  universe *u = segment_graph(width*height, num, edges, c);
  
  // post process small components
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  delete [] edges;

  image<rgb> *output = new image<rgb>(width, height);
  rgb *colors = new rgb[width*height];
/*  image<int> *output = new image<int>(width, height);*/

  // pick random colors for each component
 if (*PpmOption){
  for (int i = 0; i < width*height; i++)
    colors[i] = random_rgb();
 }
 

  for (int y = 0; y < height-1; y++) {
    int YtimesWidth = y * width;
    int temp = y;
    for (int x = 0; x < width-1; x++) {
      int comp = u->find(YtimesWidth + x);
      if (*PpmOption){
         imRef(output, x, y) = colors[comp];
      }
      SegOut[ temp] = (double) comp;
      temp += height;
    }
  }  

  delete u;

  return output;
/*  return;*/
}

#endif
