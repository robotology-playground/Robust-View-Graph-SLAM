/* Standard */

#include <algorithm>
#include <cmath>
#include <complex>
//#include <cstdint>
#include <emmintrin.h>
#include <fstream>     // library that contains file input/output functions
#include <iostream>    // std::cout
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <math.h>
#include <numeric>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stack>
#include <string.h>    // std::string     
#include <utility>     // std::pair, std::make_pair
#include <vector>

/* mex */
#include "mex.h"

/* Cholmod */
#include "cholmod.h"
#include "cholmod_core.h"
#include "cholmod_function.h"

/* SuiteSparse */
#include "SuiteSparseQR.hpp"
#include "SuiteSparseQR_definitions.h"
//#include "SuiteSparseGPU_Runtime.hpp"

/* Eigen */
#include "Eigen/Sparse"
//#include "Eigen/CholmodSupport"
//#include "Eigen/Cholesky"
#include "Eigen/Dense"
//#include "Eigen/QR"

/* Definitions */
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;
typedef Eigen::MatrixXd MatXd;
typedef Eigen::VectorXd VecXd;
typedef Eigen::Triplet<double> T;

#define pi			3.141592653589793238462643383280  /* pi */
#define NIL 		-1
#define MIN(x,y) 	(((x)<(y))?(x):(y))
#define MAX(x,y) 	(((x)>(y))?(x):(y))

