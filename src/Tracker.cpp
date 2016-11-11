/*
 * Tracker.cpp
 *
 *  Created on: Aug 10, 2016
 *      Author: tabuhashim
 */

#include "Tracker.h"
extern "C" {
#include "VLFeat.cpp"
#include "calibrated_fivepoint_helper.c"
}

//namespace Tracker {

using namespace std;
using namespace cv;



Tracker::Tracker() {
	// TODO Auto-generated constructor stub

}

Tracker::~Tracker() {
	// TODO Auto-generated destructor stub
}

vector<vector<int>> Tracker::bucketFeatures(vector<Point2f> matched, float bucket_width,float bucket_height) {
	// find max values
	float u_max = 0;
	float v_max = 0;
	for (int i=0; i<matched.size(); i++){
		if (matched[i].x>u_max) u_max=matched[i].x;
		if (matched[i].y>v_max) v_max=matched[i].y;
	}

	// allocate number of buckets needed
	int32_t bucket_cols = (int32_t)floor(u_max/bucket_width)  + 1;
	int32_t bucket_rows = (int32_t)floor(v_max/bucket_height) + 1;

	// assign matches to their buckets
	vector<int> u_ind, v_ind;
	int u, v;
	vector<vector<int>> buckets;
	buckets.resize(bucket_cols*bucket_rows);
	for (int i=0; i<matched.size(); i++){
		//u.push_back(floor(matched[i].x/bucket_width));
		//v.push_back(floor(matched[i].y/bucket_height));
		u = floor(matched[i].x/bucket_width);
		v = floor(matched[i].y/bucket_height);
		buckets[v*bucket_cols+u].push_back(i);
	}
	//u_ind = u;
	//v_ind = v;

	//sort(u_ind.begin(), u_ind.end()); // have to sort it first (unique is using the difference)
	//auto last = unique(u_ind.begin(), u_ind.end());
	//u.erase(last, u_ind.end());
	//sort(v_ind.begin(), v_ind.end()); // have to sort it first (unique is using the difference)
	//auto last = unique(v_ind.begin(), v_ind.end());
	//v_ind.erase(last, v_ind.end());
	//for (int i : v)
	//    cout << i << " ";
	//cout << "\n";

	//vector<vector<int>> buckets;
	//buckets.resize(bucket_cols*bucket_rows);

	return buckets;
}

void Tracker::calibrated_fivepoint(Eigen::MatrixXd& Evec, const vector<Point2f>& kp1, const vector<Point2f>& kp2){
	/*
	 *  notes:
	 *  1. reference:
	 *
	 *       H. Stew\'enius and C. Engels and D. Nist\'er, 2006
	 *       Recent Developments on Direct Relative Orientation,
	 *       http://vis.uky.edu/~stewe/FIVEPOINT,
	 *       http://www.vis.uky.edu/~stewe/publications/stewenius_engels_nister_5pt_isprs.pdf,
	 *
	 *       Henrik Stewenius, 2005
	 *       Grobner Basis Methods for Minimal Problems in Computer Vision
	 *       PhD Thesis, Lund University, 2005
	 *       http://www.maths.lth.se/matematiklth/personal/stewe/THESIS/
	 *
	 * 2. if this implementation is too slow, please see:
	 *
	 *       Nist\'er, D., 2004
	 *       An Efficient Solution to the Five-Point Relative Pose
	 *
	 * 3. due to varying standards of Q1 and Q2, it is very possible that you
	 *    get essential matrices which are the transpose of what your expected.
	 */

	// using Eigen
	Eigen::MatrixXd A(5,9);
	Mat1f AA = Mat1f::zeros(5,9); // using OpenCV
	Eigen::MatrixXd U, S, V; // Eigen::JacobiSVD related matrices
	Mat w, u, vt; // cv::SVD related matrices
	Eigen::MatrixXd B, C, D, E, M, X;
	Eigen::MatrixXcd cB, cC, cV, SOLS; // complex matrices due to complex eigenvectors

	// create constraint matrix A
	for(int i=0; i<5; i++){
		A.coeffRef(i,0)=kp2[i].x*kp1[i].x; AA(i,0)=A.coeffRef(i,0);
		A.coeffRef(i,1)=kp2[i].x*kp1[i].y; AA(i,1)=A.coeffRef(i,1);
		A.coeffRef(i,2)=kp2[i].x*       1; AA(i,2)=A.coeffRef(i,2);
		A.coeffRef(i,3)=kp2[i].y*kp1[i].x; AA(i,3)=A.coeffRef(i,3);
		A.coeffRef(i,4)=kp2[i].y*kp1[i].y; AA(i,4)=A.coeffRef(i,4);
		A.coeffRef(i,5)=kp2[i].y*       1; AA(i,5)=A.coeffRef(i,5);
		A.coeffRef(i,6)=       1*kp1[i].x; AA(i,6)=A.coeffRef(i,6);
		A.coeffRef(i,7)=       1*kp1[i].y; AA(i,7)=A.coeffRef(i,7);
		A.coeffRef(i,8)=       1*       1; AA(i,8)=A.coeffRef(i,8);
	} // (5x9), where 5 points are being used

	// [~,~,V] = svd(Q,0); % produces the "economy size" decomposition.
	// If Q is m-by-n with m > n, then only the first n columns of U
	// are computed and S is n-by-n.
	// For m <= n, svd(X,0) is equivalent to svd(X).
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullV|Eigen::ComputeThinU);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullV|Eigen::ComputeFullU);
	U=svd.matrixU();
	S=svd.singularValues();
	V=svd.matrixV();

	//cout << "A: " << A << endl; // 5x9
	//cout << "U: " << U << endl; // 5x5
	//cout << "S: " << S << endl; // 5x5
	cout << "V: " << V << endl; // 9x9

	SVD ssvd(AA);
	//ssvd.FULL_UV;
	ssvd.compute(AA, w, u, vt, SVD::FULL_UV);

	mysvd(A, S, V);
	cout << "V: " << V << endl;

	//cout << "AA: " << AA << endl; // 5x9
	//cout << "u: " << u << endl;
	cout << "vt: " << vt << endl;
	//cout << "w: " << w << endl;

	//Mat E = Mat::zeros(9,4,CV_32FC2);
	//Mat1f E = Mat1f::zeros(9, 4);
	E=Eigen::MatrixXd::Zero(9, 4);
	for (int i=0; i<9; i++) // 0 : 9      // E = V(:,6:9); % 10x4
		for (int j=5; j<9; j++) // 5 : 9
			E.coeffRef(i,j-5)=V.coeffRef(i,j);
	//E.at<float>(i,j-5) = V.coeffRef(i,j); // 9x4
	//E(i,j-5) = V.coeffRef(i,j); // 9x4

	//cout << "E: " << E << endl; // 9x4 (matlab says 10x4)

	//double* M = (double*)calloc(10*20,sizeof(double)); // M is 10x20
	M=Eigen::MatrixXd::Zero(10,20);
	calibrated_fivepoint_helper(E,M);

	B=M.middleCols(0,10); //M.block<9, 9>(0, 0);
	C=M.middleCols(10,10); //M.block<9, 9>(0, 10);
	//std::cout << B << std::endl;
	//std::cout << C << std::endl;
	X=B.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(C);
	M=Eigen::MatrixXd::Zero(10,10);
	M.middleRows(0,3)=-X.middleRows(0,3);
	M.middleRows(3,2)=-X.middleRows(4,2);
	M.middleRows(5,1)=-X.middleRows(7,1);
	//Y.block(0,0,3,10)=-X.block(0,0,3,10);
	//Y.block(3,0,2,10)=-X.block(4,0,2,10);
	//Y.block(5,0,1,10)=-X.block(7,0,1,10);
	M.coeffRef(6,0)=1;
	M.coeffRef(7,1)=1;
	M.coeffRef(8,3)=1;
	M.coeffRef(9,6)=1;
	//std::cout << M << std::endl;

	Eigen::EigenSolver<Eigen::MatrixXd> EIG(M);
	cV=EIG.eigenvectors(); // Eigen::MatrixXcd contains complex numbers
	//std::cout << cV << std::endl;

	cB=cV.middleRows(6,3);  // contains complex numbers
	//std::cout << cB << std::endl;

	cC=Eigen::MatrixXd::Ones(3,1)*cV.middleRows(9,1);  // contains complex numbers
	//std::cout << cC << std::endl;

	SOLS=cB.cwiseQuotient(cC); // B/C
	//std::cout << SOLS << std::endl;

	Eigen::MatrixXcd cE(SOLS.rows()+1, SOLS.cols()); // concatenate vertically
	cE << SOLS,  // contains complex numbers
			Eigen::MatrixXcd::Ones(1,10);  // contains complex numbers
	//std::cout << cEvec << std::endl;

	cE=E*cE; // contains complex numbers
	//std::cout << cEvec << std::endl;

	B=cE.cwiseAbs2(); // real only
	//std::cout << B << std::endl;

	B=B.cwiseSqrt(); // real only
	//std::cout << B << std::endl;

	B=Eigen::MatrixXd::Ones(9,1)*B.colwise().sum(); // real only
	//std::cout << B << std::endl;

	cB=cE.cwiseQuotient(B);
	//std::cout << cEvec << std::endl;

	B=cB.imag();
	C=cB.real();
	//B=B.colwise().sum();
	//std::cout << B << std::endl;
	D=Eigen::MatrixXd::Zero(cE.rows(),cE.cols());
	int k = 0;
	for (int i=0; i<B.cols();i++){
		A = B.col(i);
		if (A.isZero())
			if (Evec.cols()<1){
				D.col(k) = C.col(i);
				k+=1;
			}
	}
	Evec=D.middleCols(0,k);
	// DONE

}

void Tracker::fundamentalMatrix(Mat F, const vector<KeyPoint> kp1, const vector<KeyPoint> kp2, const Mat active) {

	// number of active p_matched
	int32_t N = kp1.size();//active.size();

	// create constraint matrix A
	Mat A(N,9,CV_32FC2);
	uchar* p;
	for(int i=0; i<A.rows; i++){
		p = A.ptr<uchar>(i);
		p[0] = kp2[i].pt.x*kp1[i].pt.x;
		p[1] = kp2[i].pt.x*kp1[i].pt.y;
		p[2] = kp2[i].pt.x;
		p[3] = kp2[i].pt.y*kp1[i].pt.x;
		p[4] = kp2[i].pt.y*kp1[i].pt.y;
		p[5] = kp2[i].pt.y;
		p[6] = kp1[i].pt.x;
		p[7] = kp1[i].pt.y;
		p[8] = 1;
	}

	// compute singular value decomposition of A
	SVD svd(A);
	Mat w, u, vt;
	//SVD::compute(A, w, u, vt);
	//R = svd.u * Mat(W) * svd.vt; //HZ 9.19
	//t = svd.u.col(2); //u3
	//Matrix U,W,V;
	//A.svd(U,W,V);

	// extract fundamental matrix from the column of V corresponding to the smallest singular value
	//F = Matrix::reshape(V.getMat(0,8,8,8),3,3);

	// enforce rank 2
	//F.svd(U,W,V);
	//W.val[2][0] = 0;
	//F = U*Matrix::diag(W)*~V;
}

/* generate image points */
void Tracker::get_aligned_point_matches(vector<point_2d>& impoints, int ncams, Mat image) {
	//vector<point_2d> impoints;
	vector<double> x, y;
	vector<int> status;

	/* get_features */
	// Detect keypoints and compute descriptors using AKAZE
	//vector<KeyPoint> kpts;
	//Mat desc;
	//Ptr<AKAZE> akaze = AKAZE::create();
	//akaze(image, noArray(), kpts, desc);

	// Use brute-force matcher to find 2-nn matches
	//We use Hamming distance, because AKAZE uses binary descriptor by default.
	//BFMatcher matcher(NORM_HAMMING);
	//vector< vector<DMatch> > nn_matches;
	//matcher.knnMatch(desc, desc, nn_matches, 2);

	// Use 2-nn matches to find correct keypoint matches
	//for(size_t i = 0; i < nn_matches.size(); i++) {
	//	DMatch first = nn_matches[i][0];
	// 	float dist1 = nn_matches[i][0].distance;
	// 	float dist2 = nn_matches[i][1].distance;
	//	if(dist1 < nn_match_ratio * dist2) {
	//		matched1.push_back(kpts1[first.queryIdx]);
	//		matched2.push_back(kpts2[first.trainIdx]);
	// 	}
	//}

	// Output results
	//Mat res;
	//drawMatches(img1, inliers1, img2, inliers2, good_matches, res);
	//imwrite("res.png", res);

	//	int npts = 100 ;
	//	for (int j=0; j<npts; j++){
	//		x.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised x points
	//		y.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised y points
	//		status.push_back( 1 ) ;
	//	}
	//	impoints.push_back(point_2d( x, y, status, 1 ));
	//	x.clear() ;
	//	y.clear() ;
	//	status.clear() ;
	//
	//	if (impoints.size()>1){
	//		/* track_features */
	//		npts = impoints[0].x.size() ;
	//		for (int i=2; i<ncams; i++){
	//			for (int j=0; j<npts; j++){
	//				x.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised x points
	//				y.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised y points
	//				status.push_back( 1 ) ;
	//			}
	//			// insert tracking code here
	//			impoints.push_back(point_2d( x, y, status, i ));
	//			x.clear() ;
	//			y.clear() ;
	//			status.clear() ;
	//		}
	//	}
	//return impoints ;
}

/* generates linearisation point */
double* Tracker::generate_linearisation_point(vector<Tracker::point_2d> p, double *Pkin) {
	int ncams, npts;
	ncams = p.size();
	npts = p[0].x.size();
	double *xs;
	//Pkin = generate_iCub_kinematics ( ncams ); /* generate forward kinematics matrices */
	//rho = generate_inverse_depth ( ncams ); /* generate point inverse depth */
	int i = 0 ;
	xs = new double [6*ncams+npts]();
	for (i=6*1; i<6*ncams; i=i+12)                 /* stereo constraints */
		xs[i] = 0.068; // Pkin goes here
	for (i=6*2; i<6*ncams; i=i+12)              /* monocular constraints */
		//xs[i] = fRand(0.03, 0.05); // Pkin goes here
		for (i=6*ncams; i<6*ncams+npts; i++)     /* inverse depth parameters */
			//xs[i] = fRand(1, 2); // rho goes here
			return xs;
}

boost::dynamic_bitset<> Tracker::remove_points_at_infinity(vector<Point2f> matched1, vector<Point2f> matched2, double thres) {
	/* matches should move at least pixel_disparity pixels
	 * the largest the pixel_disparity threshold is, the more stable the results
	 * and the estimates, but the less dense the result is and the shorter the
	 * observed range is.
	 */
	float diff1, diff2;
	boost::dynamic_bitset<> vis(matched1.size()); // all 0's by default
	for (int i=0; i<matched1.size(); i++){
		diff1 = matched2[i].x - matched1[i].x;
		diff2 = matched2[i].y - matched1[i].y;
		if ( fabs(diff1) > thres | fabs(diff2) > thres )
			vis[i] = 1;
	}
	return vis;
}

void Tracker::setFirstFrame(const Mat frame) {
	first_frame = frame.clone();
	detector->detect(first_frame, first_kp, noArray());
	cout << "Key-points : " << (int)first_kp.size() << ", " ;
	descriptor->compute(first_frame, first_kp, first_desc);
	cout << "Descriptors : " << (int)first_desc.rows << ", " << (int)first_desc.cols << endl;
}

boost::dynamic_bitset<> Tracker::verify_point_track(vector<Point2f> matched1, vector<Point2f> matched2, vector<int> tracked, vector<double> error){
	// minimum features disparity and validity
	double minerror = 1000; // very high, because is not considered at this stage
	double mindisp = 5;

	boost::dynamic_bitset<> status(matched1.size()); // all 0's by default

	//vis = remove_points_at_infinity(p1,p2,options.mindisp);
	////use = use & (vis' & status & isfinite(cr{ii}(:,1)));
	//status = vis' & isfinite(p2(:,1)) & isfinite(p2(:,2)) & p2(:,1)<options.imgsize(2) & p2(:,2)<options.imgsize(1) & p2(:,1)>0 & p2(:,2)>0 ;
	//		if nargin>3
	//			status=status&tracked;
	//		end
	//		if nargin>4
	//			status=status&error<minerror;
	//		end
	return status;
}


// Given a matrix M/A[1..m][1..n], this routine computes its singular value decomposition, M/A =
// U·W·V T. Thematrix U replaces a on output. The diagonal matrix of singular values W is output
// as a vector w[1..n]. Thematrix V (not the transpose V T ) is output as v[1..n][1..n].
void Tracker::mysvd(Eigen::MatrixXd& U, Eigen::MatrixXd& S, Eigen::MatrixXd& V) {


	int m = U.rows();
	int n = U.cols();
	//U2 = Matrix(m,m);
	//V  = Matrix(n,n);

	float* w   = (float*)malloc(n*sizeof(float));
	float* rv1 = (float*)malloc(n*sizeof(float));

	int32_t flag,i,its,j,jj,k,l,nm;
	float   anorm,c,f,g,h,s,scale,x,y,z;

	g = scale = anorm = 0.0; // Householder reduction to bidiagonal form.
	for (i=0;i<n;i++) {
		l = i+1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(U.coeffRef(k,i));
			if (scale) {
				for (k=i;k<m;k++) {
					U.coeffRef(k,i) /= scale;
					s += U.coeffRef(k,i)*U.coeffRef(k,i);
				}
				f = U.coeffRef(i,i);
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				U.coeffRef(i,i) = f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += U.coeffRef(k,i)*U.coeffRef(k,j);
					f = s/h;
					for (k=i;k<m;k++) U.coeffRef(k,j) += f*U.coeffRef(k,i);
				}
				for (k=i;k<m;k++) U.coeffRef(k,i) *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if (i<m && i!=n-1) {
			for (k=l;k<n;k++) scale += fabs(U.coeffRef(i,k));
			if (scale) {
				for (k=l;k<n;k++) {
					U.coeffRef(i,k) /= scale;
					s += U.coeffRef(i,k)*U.coeffRef(i,k);
				}
				f = U.coeffRef(i,l);
				g = -SIGN(sqrt(s),f);
				h = f*g-s;
				U.coeffRef(i,l) = f-g;
				for (k=l;k<n;k++) rv1[k] = U.coeffRef(i,k)/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += U.coeffRef(j,k)*U.coeffRef(i,k);
					for (k=l;k<n;k++) U.coeffRef(j,k) += s*rv1[k];
				}
				for (k=l;k<n;k++) U.coeffRef(i,k) *= scale;
			}
		}
		anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) { // Accumulation of right-hand transformations.
		if (i<n-1) {
			if (g) {
				for (j=l;j<n;j++) // Double division to avoid possible underflow.
					V.coeffRef(j,i)=(U.coeffRef(i,j)/U.coeffRef(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += U.coeffRef(i,k)*V.coeffRef(k,j);
					for (k=l;k<n;k++) V.coeffRef(k,j) += s*V.coeffRef(k,i);
				}
			}
			for (j=l;j<n;j++) V.coeffRef(i,j) = V.coeffRef(j,i) = 0.0;
		}
		V.coeffRef(i,i) = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i=IMIN(m,n)-1;i>=0;i--) { // Accumulation of left-hand transformations.
		l = i+1;
		g = w[i];
		for (j=l;j<n;j++) U.coeffRef(i,j) = 0.0;
		if (g) {
			g = 1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += U.coeffRef(k,i)*U.coeffRef(k,j);
				f = (s/U.coeffRef(i,i))*g;
				for (k=i;k<m;k++) U.coeffRef(k,j) += f*U.coeffRef(k,i);
			}
			for (j=i;j<m;j++) U.coeffRef(j,i) *= g;
		} else for (j=i;j<m;j++) U.coeffRef(j,i)=0.0;
		++U.coeffRef(i,i);
	}
	for (k=n-1;k>=0;k--) { // Diagonalization of the bidiagonal form: Loop over singular values,
		for (its=0;its<30;its++) { // and over allowed iterations.
			flag = 1;
			for (l=k;l>=0;l--) { // Test for splitting.
				nm = l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) { flag = 0; break; }
				if ((float)(fabs( w[nm])+anorm) == anorm) { break; }
			}
			if (flag) {
				c = 0.0; // Cancellation of rv1[l], if l > 1.
				s = 1.0;
				for (i=l;i<=k;i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y = U.coeffRef(j,nm);
						z = U.coeffRef(j,i);
						U.coeffRef(j,nm) = y*c+z*s;
						U.coeffRef(j,i)  = z*c-y*s;
					}
				}
			}
			z = w[k];
			if (l==k) { // Convergence.
				if (z<0.0) { // Singular value is made nonnegative.
					w[k] = -z;
					for (j=0;j<n;j++) V.coeffRef(j,k) = -V.coeffRef(j,k);
				}
				break;
			}
			if (its == 29)
				cerr << "ERROR in SVD: No convergence in 30 iterations" << endl;
			x = w[l]; // Shift from bottom 2-by-2 minor.
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.0; // Next QR transformation:
			for (j=l;j<=nm;j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x = V.coeffRef(jj,j);
					z = V.coeffRef(jj,i);
					V.coeffRef(jj,j) = x*c+z*s;
					V.coeffRef(jj,i) = z*c-x*s;
				}
				z = pythag(f,h);
				w[j] = z; // Rotation can be arbitrary if z = 0.
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y = U.coeffRef(jj,j);
					z = U.coeffRef(jj,i);
					U.coeffRef(jj,j) = y*c+z*s;
					U.coeffRef(jj,i) = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}

	// sort singular values and corresponding columns of u and v
	// by decreasing magnitude. Also, signs of corresponding columns are
	// flipped so as to maximize the number of positive elements.
	int32_t s2,inc=1;
	float   sw;
	float* su = (float*)malloc(m*sizeof(float));
	float* sv = (float*)malloc(n*sizeof(float));
	do { inc *= 3; inc++; } while (inc <= n);
	do {
		inc /= 3;
		for (i=inc;i<n;i++) {
			sw = w[i];
			for (k=0;k<m;k++) su[k] = U.coeffRef(k,i);
			for (k=0;k<n;k++) sv[k] = V.coeffRef(k,i);
			j = i;
			while (w[j-inc] < sw) {
				w[j] = w[j-inc];
				for (k=0;k<m;k++) U.coeffRef(k,j) = U.coeffRef(k,j-inc);
				for (k=0;k<n;k++) V.coeffRef(k,j) = V.coeffRef(k,j-inc);
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for (k=0;k<m;k++) U.coeffRef(k,j) = su[k];
			for (k=0;k<n;k++) V.coeffRef(k,j) = sv[k];
		}
	} while (inc > 1);
	for (k=0;k<n;k++) { // flip signs
		s2=0;
		for (i=0;i<m;i++) if (U.coeffRef(i,k) < 0.0) s2++;
		for (j=0;j<n;j++) if (V.coeffRef(j,k) < 0.0) s2++;
		if (s2 > (m+n)/2) {
			for (i=0;i<m;i++) U.coeffRef(i,k) = -U.coeffRef(i,k);
			for (j=0;j<n;j++) V.coeffRef(j,k) = -V.coeffRef(j,k);
		}
	}

	// create vector and copy singular values
	//W = Matrix(min(m,n),1,w);
	S=Eigen::VectorXd::Zero(IMIN(m,n));
	k=0;
	for (int32_t i=0; i<m; i++)
		S.coeffRef(i) = w[k++];

	// extract mxm submatrix U
	//U2.setMat(U.getMat(0,0,m-1,min(m-1,n-1)),0,0);

	// release temporary memory
	free(w);
	free(rv1);
	free(su);
	free(sv);
}

float Tracker::pythag(float a,float b) {
	float absa,absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0+SQR(absb/absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void Tracker::triangulate_inverse_depth(vector<double>& r, const Tracker::point_2d &p1, const Tracker::point_2d &p2, const Mat R, const Mat t){

    /* create the output matrix */
    //vector<point_3d> xf;

    double a, b, c, d, e, f, denom, s, h;
    double D, r1, r2, rim, z;
    double tx=t.at<double>(0), ty=t.at<double>(1), tz=t.at<double>(2);
    double r00=R.at<double>(0,0), r01=R.at<double>(0,1), r02=R.at<double>(0,2);
    double r10=R.at<double>(1,0), r11=R.at<double>(1,1), r12=R.at<double>(1,2);
    double r20=R.at<double>(2,0), r21=R.at<double>(2,1), r22=R.at<double>(2,2);
    double v10, v11, v12, v20, v21, v22;
    double x0, y0, x1, y1;

    for (int i=0;i<p1.x.size();i++){
    	x0 = p1.x[i];
        y0 = p1.y[i];
        x1 = p2.x[i];
        y1 = p2.y[i];

        // calculate lines
        v10 = x0;
        v11 = y0;
        v12 = 1;
        v20 = r00*x1+r01*y1+r02;
        v21 = r10*x1+r11*y1+r12;
        v22 = r20*x1+r21*y1+r22;

        // distance formulas based on Schneider pp 409-412
		  // shortest distance between the two lines
        a = v10*v10+v11*v11+v12*v12;
        b = v10*v20+v11*v21+v12*v22;
        c = v20*v20+v21*v21+v22*v22;
        d = tx*v10+ty*v11+tz*v12;
        e = tx*v20+ty*v21+tz*v22;
        f = tx*tx+ty*ty+tz*tz;
        denom = a*c-b*b;
        if (denom < std::numeric_limits<double>::epsilon()) {
            denom = 1;} // accounts for parallel lines
        s = (c*d-b*e)/denom;
        h = (b*d-a*e)/denom;

        D = f+s*(a*s-b*h-2*d)+h*(c*h-b*s+2*e);

        // compute lines length
        r1 = sqrt(a)*s;
        r2 = sqrt(b)*h;

        // negative distance
        //vis = r1 > 0;

        // inverse depth
        rim = sqrt(x0*x0+y0*y0+1);
        r.push_back(r1/rim);
    }
}

Mat Tracker::process(const Mat frame){

	/*
	// VLSIFT
	Mat temp = frame.clone();
	double* TFrames = (double*)calloc ( 4 * 10000, sizeof(double) );
	double* QFrames = (double*)calloc ( 4 * 10000, sizeof(double) );
	uint8_t* TDescr  = (uint8_t*)calloc ( 128 * 10000, sizeof(uint8_t) );
	uint8_t* QDescr  = (uint8_t*)calloc ( 128 * 10000, sizeof(uint8_t) );
	int Tnframes = 0;
	int Qnframes = 0;
	VLSIFT(temp, TDescr2, TFrames, &Tnframes);
	VLSIFT(temp, QDescr2, QFrames, &Qnframes);
	cout << "VLSIFT : " << Tnframes << " " << Qnframes << endl;
	TFrames=(double*)realloc(TFrames,4*sizeof(double)*Tnframes); // = Y X Scale Angle
	TDescr=(uint8_t*)realloc(TDescr,128*sizeof(uint8_t)*Tnframes);
	 */

	//-- Detect and compute
	Mat second_frame = frame.clone();
	vector<KeyPoint> kp;
	Mat desc;
	detector->detect(second_frame, kp, noArray());
    cout << "Key-points : " << (int)kp.size() << ", " ;
	descriptor->compute(second_frame, kp, desc);
    cout << "Descriptors : " << (int)desc.rows << ", " << (int)desc.cols << endl;

	/*
	//-- Moves from vector<KeyPoint> to double*
	double* TFrames = (double*)calloc(2*kp.size(),sizeof(double));
	double* QFrames = (double*)calloc(2*first_kp.size(),sizeof(double));
	for(int i=0; i<first_desc.rows; i++){
		QFrames[2*i+0] = first_kp[i].pt.x;
		QFrames[2*i+1] = first_kp[i].pt.y;
	}
	for(int i=0; i<desc.rows; i++){
		TFrames[2*i+0] = kp[i].pt.x;
		TFrames[2*i+1] = kp[i].pt.y;
	}
	 */

	//-- Moves from Mat to double*
	uint8_t* TDescr = (uint8_t*)calloc(128*desc.rows,sizeof(uint8_t));
	uint8_t* QDescr = (uint8_t*)calloc(128*first_desc.rows,sizeof(uint8_t));
	uchar* p;
	for(int i=0; i<first_desc.rows; i++){
		p = first_desc.ptr<uchar>(i);
		for(int j=0; j<first_desc.cols; j++)
			QDescr[first_desc.cols*i+j] = p[j];
	}
	for(int i=0; i<desc.rows; i++){
		p = desc.ptr<uchar>(i);
		for(int j=0; j<desc.cols; j++)
			TDescr[desc.cols*i+j] = p[j];
	}

	//-- VLMATCH
	int Tnframes=(int)kp.size();
	int Qnframes=(int)first_kp.size();
	int matchesFound=0;
	double* MATCHES=(double*)calloc(10000,sizeof(double));
	VLMATCH(TDescr, QDescr, Tnframes, Qnframes, .85, &matchesFound, MATCHES);
	MATCHES=(double*)realloc(MATCHES,sizeof(double)*matchesFound*2);
	cout << "VLFEAT Matches : " << matchesFound << endl;

	//-- Opencv matching
	vector<DMatch> matches;
	Mat mask;
	matcher->match(first_desc, desc, matches, noArray());
	cout << "OpenCV Matches : " << (int)matches.size() << endl;

	double max_dist = 0;
	double min_dist = 100;
	double mil_dist = 5; // distance multiplier

	//-- Quick calculation of max and min distances between keypoints
	// http://docs.opencv.org/3.1.0/d5/d6f/tutorial_feature_flann_matcher.html#gsc.tab=0
	for(int i = 0;i<matches.size();i++){
		double dist = matches[i].distance;
		if(dist<min_dist) min_dist=dist;
		if(dist>max_dist) max_dist=dist;
	}
	//cout<<"-- Max dist : " <<max_dist<<"\n";
	//cout<<"-- Min dist : " <<min_dist<<"\n";

	//-- Get only "good" matches (i.e. whose distance is less than 2*min_dist,
	//-- or a small arbitary value ( 0.02 ) in the event that min_dist is very
	//-- small)
	//-- PS.- radiusMatch can also be used here.
	vector<DMatch> good_matches;
	for(int i=0;i<matches.size();i++)
		if(matches[i].distance<=max(mil_dist*min_dist,0.02))
			good_matches.push_back(matches[i]);

	//-- Collect the matching key-points
	vector<Point2f> matched1, matched2;
	//for(int i = 0; i < matches.size(); i++){ // using OpenCV matches
	//	matched1.push_back(first_kp[matches[i].queryIdx].pt);
	///	matched2.push_back(      kp[matches[i].trainIdx].pt);
	//}
	for(int i=0;i<matchesFound;i++){ // using VLMATCH
		matched1.push_back(first_kp[MATCHES[2*i+0]].pt);
		matched2.push_back(      kp[MATCHES[2*i+1]].pt);
	}

	//-- Remove small tracks
	boost::dynamic_bitset<> vis;
	vis = remove_points_at_infinity(matched1, matched2, 2);
	vector<Point2f> vismatched1, vismatched2;
	for(int i = 0; i < matches.size(); i++)
		if(vis[i]){
			vismatched1.push_back(matched1[i]);
			vismatched2.push_back(matched2[i]);
		}
	cout << "Visible: " << (int)vis.count() << endl;

	// drawing the results
	namedWindow("matches", 1);
	Mat img_matches;
	drawMatches(first_frame, first_kp, frame, kp, good_matches, img_matches);
	imshow("matches", img_matches);
	waitKey(0);

	if (vis.count()>100){
		//-- http://nghiaho.com/?p=1675
		//-- https://avisingh599.github.io/vision/monocular-vo/
		Mat K1 = (Mat_<double>(3,3) << 235.1162, 0, 155.6680, 0, 235.7933, 122.0000, 0, 0, 1);
		Mat K2 = (Mat_<double>(3,3) << 234.2173, 0, 149.1942, 0, 234.8432, 124.5134, 0, 0, 1);
		Mat k1 = (Mat_<double>(1,5) << -0.4323, 0.1952, -0.0003,  0.0018, 0);
		Mat k2 = (Mat_<double>(1,5) << -0.4309, 0.1880,  0.0007, -0.0008, 0);

		//vector<vector<int>> buckets;
		//buckets = bucketFeatures(vismatched2, 10, 10);
		//for (int i=0; i<buckets.size(); i++){
		//	for (int j=0; j<buckets[i].size(); j++)
		//		cout << buckets[i][j] << " ";
		//	if (buckets[i].size()>1)
		//		cout << endl;
		//}

		//-- Normalize
		vector<Point2f> vismatched1_, vismatched2_;
		undistortPoints(vismatched1, vismatched1_, K1, k1, noArray(), K1);
		undistortPoints(vismatched2, vismatched2_, K2, k2, noArray(), K2);

		//-- Essential matrix with RANSAC
		Mat E, R, t, mask;
		E = findEssentialMat(vismatched2, vismatched1, 1.0, Point2d(0,0), RANSAC, 0.999, 0.0001, mask);
		//correctMatches(E, vismatched1, vismatched2, vismatched1, vismatched2);

		//Mat F = findFundamentalMat(vismatched2_, vismatched1_, FM_RANSAC, 0.1, 0.99, mask);
		//correctMatches(F, vismatched1_, vismatched2_, vismatched1_, vismatched2_);
		//E = K2.t()*F*K1;

		//SVD svd(E);
		//Matx33d W(0,-1,0,   //HZ 9.13
		//      1,0,0,
		//      0,0,1);
		//Matx33d Winv(0,1,0,
		//     -1,0,0,
		//     0,0,1);
		//R = svd.u * Mat(W) * svd.vt; //HZ 9.19
		//t = svd.u.col(2); //u3

		//-- Pose recovery
		undistortPoints(vismatched1, vismatched1_, K1, k1);
		undistortPoints(vismatched2, vismatched2_, K2, k2);
		recoverPose(E, vismatched2_, vismatched1_, R, t, 1.0, Point2d(0,0), noArray());

		Mat a;
		Rodrigues(R,a);
		transpose(a,a);
		transpose(t,t);
		cout << t << " " << a*180/M_PI << endl ;

		//vector<Point2f> kp1, kp2;
		//for (int i=0; i<5; i++){
		//	kp1.push_back(vismatched1_[i]);
		//	kp2.push_back(vismatched2_[i]);
		//}

		//kp1[0].x=0.7060;kp1[1].x=0.0971;kp1[2].x=0.9502;kp1[3].x=0.7655;kp1[4].x=0.4456;
		//kp1[0].y=0.0318;kp1[1].y=0.8235;kp1[2].y=0.0344;kp1[3].y=0.7952;kp1[4].y=0.6463;
		//kp2[0].x=0.2769;kp2[1].x=0.6948;kp2[2].x=0.4387;kp2[3].x=0.1869;kp2[4].x=0.7094;
		//kp2[0].y=0.0462;kp2[1].y=0.3171;kp2[2].y=0.3816;kp2[3].y=0.4898;kp2[4].y=0.7547;

		Eigen::MatrixXd Evec, E2(3,3);
		calibrated_fivepoint(Evec, vismatched1_, vismatched2_);
		//cout << Evec << endl;
		//cout << endl;
		for (int j=0; j<Evec.cols(); j++){
			for (int i=0; i<Evec.rows(); i++){
				int l = floor(i/3);
				int k = i - (l-1)*3 - 3;
                //cout << k << " " << l << endl;
				//E.at<double>(k,l) = Evec.coeffRef(i,j)/Evec.coeffRef(8,j); // normalise to force rank 2
				E2.coeffRef(k,l) = Evec.coeffRef(i,j)/Evec.coeffRef(8,j);
			}
			//cout << E << endl;
			//cout << endl;

			//Mat U, W, Vt; // cv::SVD related matrices
			//SVD svd(E);
			//svd.compute(E, W, U, Vt, SVD::FULL_UV);

			Eigen::JacobiSVD<Eigen::MatrixXd> svd1(E2,Eigen::ComputeFullV|Eigen::ComputeFullU);
			auto U = svd1.matrixU();
			auto S = svd1.singularValues();
			auto V = svd1.matrixV();

			// Generate motion hypothesis
			// SVD decomposition of the essential matrix
			// re-enforce rank 2 constraint on essential matrix
			//Eigen::DiagonalMatrix<double,3> D;
			S.coeffRef(2) = 0 ;
			auto D = S.asDiagonal() ;
			E2 = U*D*V.transpose() ; //E = U*diag([D(1,1) D(2,2) 0])*V';
			// SVD decomposition of the essential matrix
			Eigen::JacobiSVD<Eigen::MatrixXd> svd2(E2,Eigen::ComputeFullV|Eigen::ComputeFullU) ;
			U = svd2.matrixU() ;
			V = svd2.matrixV() ;
			Eigen::MatrixXd W = U*V.transpose() ;
			if ( W.determinant()<0 ) V = -V ;//(det(U*V')<0)
			W = Eigen::MatrixXd::Zero(3,3) ;
			W.coeffRef(0,1) = -1 ;
			W.coeffRef(1,0) = +1 ;
			W.coeffRef(2,2) = +1 ;
			//W = [0 -1 0; +1 0 0 ; 0 0 +1]; % Harley's book
			//Z = [0 +1 0; -1 0 0 ; 0 0  0]; % Harley's book
			Eigen::MatrixXd R1 = U*W*V.transpose() ;
			Eigen::MatrixXd R2 = U*W.transpose()*V.transpose() ;
			Eigen::MatrixXd t1 = U.middleCols(2,1) ;
			Eigen::MatrixXd t2 = -t1 ;
			//% assure determinant to be positive
			if ( R1.determinant() < 0 ) R1 = -R1 ;
			if ( R2.determinant() < 0 ) R2 = -R2 ;

			//function [P, inlier, err, sol] = resolve_motion_ambiguity(R1, R2, t1, t2, p1, p2, pixtol)
			//% case 1: {R1,t1}
			//xc = [t1; R2w(R1)'];
			//r1 = test_triangulate_inverse_depth(p1, p2, xc);
			//%xc = transform_to_relative_w(zeros(6,1), xc);
			//%xc = [-R1'*t1; R2w(R1')'];
			//%r2 = test_triangulate_inverse_depth(p2, p1, xc);
			//posDepth = sum(r1 > 0);
			//P = [R1, t1]; sol = 1;
			//% case 2: {R1,t2}
			//xc = [t2; R2w(R1)'];
			//r1 = test_triangulate_inverse_depth(p1, p2, xc);
			//%xc = transform_to_relative_w(zeros(6,1), xc);
			//%xc = [-R1'*t2; R2w(R1')'];
			//%r2 = test_triangulate_inverse_depth(p2, p1, xc);
			//if sum(r1 > 0) > posDepth
			//		P = [R1, t2]; sol = 2;
			//posDepth = sum(r1 > 0);
			//end
			//% case 3: {R2,t1}
			//xc = [t1; R2w(R2)'];
			//r1 = test_triangulate_inverse_depth(p1, p2, xc);
			//%xc = transform_to_relative_w(zeros(6,1), xc);
			//%xc = [-R2'*t1; R2w(R2')'];
			//%r2 = test_triangulate_inverse_depth(p2, p1, xc);
			//if sum(r1 > 0) > posDepth
			//		P = [R2, t1]; sol = 2;
			//posDepth = sum(r1 > 0);
			//end
			//% case 4: {R2,t2}
			//xc = [t2; R2w(R2)'];
			//r1 = test_triangulate_inverse_depth(p1, p2, xc);
			//%xc = transform_to_relative_w(zeros(6,1), xc);
			//%xc = [-R2'*t2; R2w(R2')'];
			//%r2 = test_triangulate_inverse_depth(p2, p1, xc);
			//if sum(r1 > 0) > posDepth
			//		P = [R2, t2]; sol = 2;
			//posDepth = sum(r1 > 0);
			//end
			//if posDepth
			//R = P(1:3,1:3); t = P(1:3,4);
			//% minimum depth
			//xc = [t; R2w(R)'];
			//r1 = test_triangulate_inverse_depth(p1, p2, xc);
			//%xc = transform_to_relative_w(zeros(6,1), xc);
			//%xc = [-R'*t; R2w(R')'];
			//%r2 = test_triangulate_inverse_depth(p2, p1, xc);
			//% error
			//z1 = get_scan_from_range(p1,1./r1);
			//z2 = z1;
			//z2(1,:) = z1(1,:) - t(1);
			//z2(2,:) = z1(2,:) - t(2);
			//z2(3,:) = z1(3,:) - t(3);
			//z2 = pflat(R'*z2);
			//		err = sqrt(sum((p2(1:2,:)-z2(1:2,:)).^2));
			//%err = sqrt(sum((p1(1:2,:)-z1(1:2,:)).^2) + sum((p2(1:2,:)-z2(1:2,:)).^2));
			//% inliers
			//inlier = (err < pixtol) & (r1 > 0);%&(r2 > 0);
			//end
		}

	}

	Mat res;
	return res;

}

//} /* namespace Tracker */
