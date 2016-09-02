/*
 * Tracker.h
 *
 *  Created on: Aug 10, 2016
 *      Author: tabuhashim
 */

#ifndef TRACKER_H_
#define TRACKER_H_

#include <eigens.h>
#include <opencvs.h>
#include <boost/dynamic_bitset.hpp>

#define SWAP(a,b) {temp=a;a=b;b=temp;}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
static int32_t iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

//namespace Tracker {
class Tracker {
// http://docs.opencv.org/3.0-beta/doc/tutorials/features2d/akaze_tracking/akaze_tracking.html
// http://docs.opencv.org/trunk/d8/d30/classcv_1_1AKAZE.html#aaeec869ae038190ffac3f38174058b25
public:
	Tracker();
	virtual ~Tracker();
    Tracker(cv::Ptr<cv::Feature2D> _detector, cv::Ptr<cv::Feature2D> _descriptor, cv::Ptr<cv::DescriptorMatcher> _matcher) :
		detector(_detector), descriptor(_descriptor), matcher(_matcher) {}
	void setFirstFrame(const cv::Mat frame);
	cv::Mat process(const cv::Mat frame);
	cv::Ptr<cv::Feature2D> getDetector() { return detector; }
	struct point_2d { /* image points structure */
		const std::vector<double> x ; /* x image coordinate */
		const std::vector<double> y ; /* y image coordinate */
		const std::vector<int> status ; /* status */
		int camera ; // camera id
		int npts ; // number of tracks in this image
		point_2d ( ) ; // Do not use this default constructor (structure will be empty)
		// Member initialization in a constructor
		point_2d ( std::vector<double>& a, std::vector<double>& b,
				std::vector<int>& c, int d) : x(a), y(b), status(c), camera(d) {
			npts=x.size(); }
	};
protected:
	cv::Ptr<cv::Feature2D> detector;
    cv::Ptr<cv::Feature2D> descriptor;
	cv::Ptr<cv::DescriptorMatcher> matcher;
	cv::Mat first_frame, first_desc;
	std::vector<cv::KeyPoint> first_kp;
private:
	std::vector<std::vector<int>> bucketFeatures(std::vector<cv::Point2f> matched, float bucket_width,float bucket_height);
	void calibrated_fivepoint(Eigen::MatrixXd& Evec, const std::vector<cv::Point2f>& kp1, const std::vector<cv::Point2f>& kp2);
	void fundamentalMatrix (cv::Mat F, const std::vector<cv::KeyPoint> kp1, const std::vector<cv::KeyPoint> kp2, const cv::Mat active);
	void get_aligned_point_matches(std::vector<point_2d>& impoints, int ncams, cv::Mat image);
	double* generate_linearisation_point(std::vector<point_2d> p, double *Pkin);
	boost::dynamic_bitset<> remove_points_at_infinity(std::vector<cv::Point2f> matched1, std::vector<cv::Point2f> matched2, double thres);
	boost::dynamic_bitset<> verify_point_track(std::vector<cv::Point2f> matched1, std::vector<cv::Point2f> matched2, std::vector<int> tracked, std::vector<double> error);

	float pythag(float a,float b);
	void mysvd(Eigen::MatrixXd& U, Eigen::MatrixXd& S, Eigen::MatrixXd& V);
	void triangulate_inverse_depth(std::vector<double>& r, const  Tracker::point_2d &p1, const Tracker::point_2d &p2, const cv::Mat R, const cv::Mat t);
};
//} /* namespace Tracker */

#endif /* TRACKER_H_ */
