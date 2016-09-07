#include <cpps.h>
#include <yarps.h>
#include <icubs.h>
//#include "Image.h"
//#include "ImageViewer.h"
#include "Tracker.h"
#include "PwgOptimiser.h"
#include "GraphOptimiser.h"
//#include "5point.cpp"
//#include "Rpoly.cpp"

# define	M_PI	3.14159265358979323846  /* pi */
# define    vgSLAM_KAZE    0 //det
# define    vgSLAM_FAST    1 //det
# define    vgSLAM_SIFT    2 //det & desc
# define    vgSLAM_GFTT    3 //det
# define    vgSLAM_SURF    4 //det & desc
# define    vgSLAM_BRIEF   5 //desc
# define    vgSLAM_ORB     6 //det & desc
# define    vgSLAM_BRISK   7 //det & desc
# define    vgSLAM_FREAK   8 //desc
# define    vgSLAM_FLANN   9 //matcher
# define    vgSLAM_BRUTEFORCEL2   10 //matcher
# define    vgSLAM_BRUTEFORCEL1 11 //matcher
# define    vgSLAM_BRUTEFORCEHAMMING   12//matcher
# define    vgSLAM_AKAZE    13 //det

using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;
using namespace iCub::iKin;
using namespace std;
using namespace cv;

const int num_ransac_itr = 300;
const double akaze_thresh = 1e-4; // AKAZE detection threshold set to locate about 1000 keypoints
const double ransac_thresh = 2.5f; // RANSAC inlier threshold

/*
 *you may use Valgrind to test for any memory leaks
 *Valgrind is available for download at http://valgrind.org/
 *after compiling and installation, run
 *          valgrind ./PwgOptimiser 2 2
 */

/* generates random double in the range ( fMin:fMax ) */
double fRand(double fMin, double fMax) {
	double f = (double) std::rand() / RAND_MAX ;
	return fMin + f * (fMax - fMin) ;
}

/* generates bundle constraints */
void generate_constraints_image_inverse_depth_Mviews(PwgOptimiser *Object, std::vector<Tracker::Tracker::point_2d> p){
	int ncams, npts ;
	ncams = p.size() ;
	npts = p[0].x.size() ;
	int cam, kpt, sw=0;
	std::vector<double> p1(2), z(2), R(4,0.0) ;
	Eigen::MatrixXd yz = Eigen::MatrixXd::Zero(7,1) ;
	Eigen::VectorXd Yz = Eigen::MatrixXd::Zero(7,7) ;
	for (int i=0; i<ncams; i++)
		for (int j=0; j<npts; j++)
			if (p[i].status[j]==1) {
				cam = i+1; // this should take values in the range 1:ncams
				kpt = j+1; // this should take values in the range 1:npts
				p1[0] = p[0].x[j] ; p1[1] = p[0].y[j] ;
				z[0] = p[i].x[j] ; z[1] = p[i].y[j] ;
				R[0] = 1 ; R[3] = 1 ;
				Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz, sw) ;
			}
}

/* run the vision process */
//void process (int ncams) {

	/* generate point tracks
	 * here goes:
	 *      1- Features extraction
	 *      2- Tracking or matching
	 *      3- Visibility analysis
	 *      4- Temporary/Real-time VO/kinematics solution ?
	 */
//	std::vector<Tracker::Tracker::point_2d> p ;
	//generate_point_tracks ( p, ncams ) ; // thread safe
//	ncams = p.size() ;
//	int npts = p[0].x.size() ;
	//for (int i=0; i<ncams; i++)
	//    std::cout << p[i].x[0] << " " << p[i].y[0] << " " << p[i].status[0] << std::endl ;

	/* initialise a linearisation point
	 * here goes:
	 *      1- kinematics
	 *      2- triangulation
	 *      3- Possibly output a temporary solution before the batch one?
	 *          for this one, generate_point_tracks and generate_linearisation_point
	 *          need to be merged ?????
	 *          Then, generate_point_tracks should be integrated into
	 *          generate_linearisation_point ?????
	 *          Then xs is accomulated incrementally rather than in one step
	 *          OR, better leave generate_point_tracks do that stuff alone ?????
	 */
//	double *xs, *Pkin ;
//	xs = Tracker::generate_linearisation_point( p, Pkin ) ;
	//for (int i=0; i<6*ncams+npts; i++) std::cout << xs[i] << std::endl ;

	/* initialise a PwgOptimiser object */
//	PwgOptimiser *Object ; // pointer initialisation
//	Object = new PwgOptimiser ( ncams, npts ) ; // pointer initialisation

	/* generate constraints */
//	generate_constraints_image_inverse_depth_Mviews( Object, p ) ;

	/* optimise constraints information */
//	Object->optimise_constraints_image_inverse_depth_Mviews( xs ) ;

	/* free memory */
//	delete[] xs ;
//	delete[] Pkin ;
//	delete Object ; // delete class pointer
//}
int assignMethod(string str){
    if (str.compare("KAZE")==0)
        return vgSLAM_KAZE;
    if (str.compare("FAST")==0)
        return vgSLAM_FAST;
    if (str.compare("SIFT")==0)
        return vgSLAM_SIFT;
    if (str.compare("GFTT")==0)
        return vgSLAM_GFTT;
    if (str.compare("SURF")==0)
        return vgSLAM_SURF;
    if (str.compare("BRIEF")==0)
        return vgSLAM_BRIEF;
    if (str.compare("ORB")==0)
        return vgSLAM_ORB;
    if (str.compare("BRISK")==0)
        return vgSLAM_BRISK;
    if (str.compare("FREAK")==0)
        return vgSLAM_FREAK;
    if (str.compare("FLANN")==0)
        return vgSLAM_FLANN;
    if (str.compare("BRUTEFORCEL1")==0)
        return vgSLAM_BRUTEFORCEL1;
    if (str.compare("BRUTEFORCEL2")==0)
        return vgSLAM_BRUTEFORCEL2;
    if (str.compare("BRUTEFORCEHAMMING")==0)
        return vgSLAM_BRUTEFORCEHAMMING;
    if (str.compare("AKAZE")==0)
        return vgSLAM_AKAZE;
    return -1;
}
bool checkDetector(string str){
    if(str.compare("KAZE")!=0 && str.compare("AKAZE")!=0  && str.compare("FAST")!=0 && str.compare("SIFT")!=0
            && str.compare("GFTT")!=0 && str.compare("SURF")!=0
            && str.compare("ORB")!=0 && str.compare("BRISK")!=0){
        cout<<"Wrong detector argument, available options are: KAZE, FAST, SIFT, GFTT, SURF, ORB, BRISK, AKAZE"<<endl;
        return false;
    }
    else
        return true;
}
bool checkDescriptor(string str){
    if(str.compare("SIFT")!=0 && str.compare("SURF")!=0 && str.compare("BRIEF")!=0
            && str.compare("ORB")!=0 && str.compare("BRISK")!=0 && str.compare("FREAK")!=0){
        cout<<"Wrong descriptor argument, available options are: SIFT, SURF, BRIEF, ORB, BRISK, FREAK"<<endl;
        return false;
    }
    else
        return true;
}
bool checkMatcher(string str){
    if(str.compare("FLANN")!=0 && str.compare("BRUTEFORCEL2")!=0
            && str.compare("BRUTEFORCEL1")!=0 && str.compare("BRUTEFORCEHAMMING")!=0){
        cout<<"Wrong matcher argument, available options are: FLANN, BRUTEFORCEL2, BRUTEFORCEL1, BRUTEFORCEHAMMING"<<endl;
        return false;
    }
    else
        return true;
}

/* main code */
int main (int argc, char** argv) {
	/* we need at least one input, ncams */
	// argv[0] is the program name
    // argv[1:n] are the program input arguments
    int detID=vgSLAM_KAZE, descID=vgSLAM_SIFT, matchID=vgSLAM_FLANN; // default

	if (argc<2) {
        std::cerr << "Usage: ./vgSLAM (int)ncams int(versbose)" << std::endl;
		return 1;
	}
	int ncams = std::atoi ( argv[1] );
    int VERBOSE = 0;
    ResourceFinder rf,rfdet,rfdesc,rfmatch;
    if(!rf.setDefaultConfigFile("../../conf/vgSLAM.ini"))
        return -1;
    rf.configure(argc, argv);
    rfdesc.configure(argc, argv);
    rfmatch.configure(argc, argv);

    if (argc > 2){
        VERBOSE = std::atoi ( argv[2] );
    }

    if(checkDetector(rf.find("Detector").asString()))
        detID=assignMethod(rf.find("Detector").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default detector:KAZE"<<endl;
    if(checkDescriptor(rf.find("Descriptor").asString()))
        descID=assignMethod(rf.find("Descriptor").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default descriptor:SIFT"<<endl;
    if(checkMatcher(rf.find("Matcher").asString()))
        matchID=assignMethod(rf.find("Matcher").asString());
    else
        cout<<"Error in vgSLAM.ini, setting default matcher:FLANN"<<endl;
	//if (VERBOSE==2){
	//	cvNamedWindow( "img_L", CV_WINDOW_AUTOSIZE );
	//	cvNamedWindow( "img_R", CV_WINDOW_AUTOSIZE );
	//}

	int i=1;
	bool res=false;
	bool first=true;

	// YARP
	Network yarp; // 64 bytes still reachable
	BufferedPort<ImageOf<PixelRgb> > image1_port, image2_port;
	int image1_start=0, image2_start=0;
	image1_port.open("/vgSLAM/cam/left");
	image2_port.open("/vgSLAM/cam/right");
	yarp.connect("/icub/cam/left", image1_port.getName());
	yarp.connect("/icub/cam/right",image2_port.getName());

	// initialise tracker
	std::vector<Tracker::point_2d> p ;
    // detector switch
    Ptr<Feature2D> detector, descriptor;
    switch (detID) {
    case vgSLAM_KAZE : {
        rfdet.setDefaultConfigFile("../../conf/KAZE.ini");
        rfdet.configure(argc, argv);
        detector=KAZE::create(rfdet.find("extended").asBool(),rfdet.find("upright").asBool(),
                                    (float) rfdet.find("threshold").asDouble(),
                                    rfdet.find("nOctaves").asInt(),rfdet.find("nOctaveLayers").asInt());
        //TODO MANCA LO Switch per l'ultimo argomento.
        break;
    }
    case vgSLAM_FAST : {
        rfdet.setDefaultConfigFile("../../conf/FAST.ini");
        rfdet.configure(argc, argv);
        detector=FastFeatureDetector::create(rfdet.find("threshold").asInt(),rfdet.find("nonmaxSuppression").asBool());
        //TODO manca lo switch dell'ultimo argomento
        break;
    }
    case vgSLAM_SIFT : {
        rfdet.setDefaultConfigFile("../../conf/SIFT.ini");
        rfdet.configure(argc, argv);
        detector=xfeatures2d::SIFT::create(rfdet.find("nfeatures").asInt(),rfdet.find("nOctaveLayers").asInt(),
                                                              rfdet.find("contrastThreshold").asDouble(),
                                                              rfdet.find("edgeThreshold").asDouble(),rfdet.find("sigma").asDouble());
        break;
    }
    case vgSLAM_GFTT: {
        rfdet.setDefaultConfigFile("../../conf/GFTT.ini");
        rfdet.configure(argc, argv);
        detector=GFTTDetector::create(rfdet.find("maxCorners").asInt(),rfdet.find("qualityLevel").asDouble(),
                                                    rfdet.find("minDistance").asDouble(),rfdet.find("blockSize").asInt(),
                                                    rfdet.find("useHarrisDetector").asBool(),rfdet.find("k").asDouble());


        break;
    }
    case vgSLAM_SURF : {
        rfdet.setDefaultConfigFile("../../conf/SURF.ini");
        rfdet.configure(argc, argv);
        detector=xfeatures2d::SURF::create(rfdet.find("hessianThreshold").asDouble(),rfdet.find("nOctaves").asInt(),
                                                              rfdet.find("inOctaveLayers").asInt(),rfdet.find("extended").asBool(),
                                                              rfdet.find("upright").asBool());
        break;
    }
    case vgSLAM_ORB : {
        rfdet.setDefaultConfigFile("../../conf/ORB.ini");
        rfdet.configure(argc, argv);
        detector=ORB::create(rfdet.find("nfeatures").asInt(),rfdet.find("scaleFactor").asDouble(),
                                 rfdet.find("nlevels").asInt(),rfdet.find("edgeThreshold").asInt(),
                                 rfdet.find("firstLevel").asInt(),rfdet.find("WTA_K").asInt(),
                                 ORB::HARRIS_SCORE,rfdet.find("patchSize").asInt(),
                                 rfdet.find("fastThreshold").asInt());
        //TODO Argomento
        break;
    }
    case vgSLAM_BRISK : {
        rfdet.setDefaultConfigFile("../../conf/BRISK.ini");
        rfdet.configure(argc, argv);
        detector=BRISK::create(rfdet.find("thresh").asInt(),rfdet.find("octaves").asInt(),rfdet.find("patternScale").asDouble());
        break;
    }
    case vgSLAM_AKAZE : {
        rfdet.setDefaultConfigFile("../../conf/AKAZE.ini");
        rfdet.configure(argc, argv);
        detector=AKAZE::create(AKAZE::DESCRIPTOR_MLDB,rfdet.find("descriptor_size").asInt(),rfdet.find("descriptor_channels").asInt(),
                               rfdet.find("threshold").asDouble(),rfdet.find("nOctaves").asInt(),rfdet.find("nOctaveLayers").asInt(),
                               KAZE::DIFF_PM_G2);
        //TODO primo e ultimo argomento
        break;
    }
    default:
        break;
    }
    // descriptor switch
    switch (descID) {
    case vgSLAM_SIFT : {
        descriptor=xfeatures2d::SIFT::create();
        break;
    }
    case vgSLAM_SURF : {
        descriptor=xfeatures2d::SURF::create();
        break;
    }
    case vgSLAM_BRIEF : {
        descriptor=xfeatures2d::BriefDescriptorExtractor::create();
        break;
    }
    case vgSLAM_ORB : {
        descriptor=ORB::create();
        break;
    }
    case vgSLAM_BRISK : {
        descriptor=BRISK::create();
        break;
    }
    case vgSLAM_FREAK : {
        descriptor=xfeatures2d::FREAK::create();
        break;
    }
    default:
        break;
    }
    // matcher switch
    Ptr<DescriptorMatcher> matcher;
    switch (matchID) {
    case vgSLAM_FLANN:{
        matcher=DescriptorMatcher::create("FlannBased");
        break;
    }
    case vgSLAM_BRUTEFORCEL2:{
        matcher=DescriptorMatcher::create("BruteForce");
        break;
    }
    case vgSLAM_BRUTEFORCEL1:{
        matcher=DescriptorMatcher::create("BruteForce-L1");
        break;
    }
    case vgSLAM_BRUTEFORCEHAMMING:{
        matcher=DescriptorMatcher::create("BruteForce-Hamming");
        break;
    }
    default:
        break;
    }
    //detector->setThreshold(akaze_thresh);//TODO how we manage different detector and different threshold and params

	//int nimages = 2*5063;
	while(i<ncams){ //i<nimages

		ImageOf<PixelRgb> *image1_yarp = image1_port.read();
		ImageOf<PixelRgb> *image2_yarp = image2_port.read();
		Stamp s1,s2;

		if(image1_port.getEnvelope(s1) && image2_port.getEnvelope(s2)){
			if(i==1){
				image1_start = s1.getCount();
				image2_start = s2.getCount();
			}
			if(abs(s1.getCount()-image1_start)>2 || abs(s2.getCount()-image2_start)>2
					|| fabs((s1.getTime())-(s2.getTime()))>0.03){//0.03 is the half delta t
				image1_start = s1.getCount();
				image2_start = s2.getCount();
				continue;
			}
			if (image1_yarp!=NULL && image2_yarp!=NULL){
				std::cout << "[" << i << "," << i+1 << "]" << std::endl ;

				/* the images */
				Mat image1_cv = cvarrToMat(static_cast<IplImage*>(image1_yarp->getIplImage()));
				cvtColor(image1_cv, image1_cv, CV_RGB2BGR);
				cvtColor(image1_cv, image1_cv, COLOR_BGR2GRAY);
				Mat image2_cv = cvarrToMat(static_cast<IplImage*>(image2_yarp->getIplImage()));
				cvtColor(image2_cv, image2_cv, CV_RGB2BGR);
				cvtColor(image2_cv, image2_cv, COLOR_BGR2GRAY);

				/* the tracker */
                Tracker tracker(detector, descriptor, matcher); // a tracker for each key-frame
				//if(first){
                tracker.setFirstFrame(image1_cv);
				//	first = false;
				//}
				//else
                //	tracker.process(image1_cv);
                tracker.process(image2_cv);

				//get_aligned_point_matches (p, ncams, image1_cv) ;// thread safe
				//get_aligned_point_matches (p, ncams, image2_cv) ;// thread safe

				//ncams = p.size() ;
				//int npts = p[0].x.size() ;

				if (VERBOSE==2){
					//cvShowImage( "img_L", image1_cv ); // Large memory leak
					//cvShowImage( "img_R", image2_cv ); // Large memory leak
					//cvWaitKey(1);                     // Large memory leak
					//cvDestroyWindow( "img_L" );
					//cvDestroyWindow( "img_R" );
				}

				//cvReleaseImage( &cvImageL );
				//cvReleaseImage( &cvImageR );
			}

			i = i + 2;
		}
	}

	// run the process
	//process( ncams ) ;

	return 0 ;
}



/*
		int num_pts = 100;
		vector <Mat> Evec; // essential matrix
		vector <Mat> Pvec; // 3x4 projection matrix
	    vector <int> inliers;
		bool ret;
		double *pts1 = new double [2*matched1.size()];
		double *pts2 = new double [2*matched2.size()];
		for(int i=0; i<matched1.size(); i++){
			pts1[2*i+0] = matched1[i].x;
			pts1[2*i+1] = matched1[i].y;
			pts2[2*i+0] = matched2[i].x;
			pts2[2*i+1] = matched2[i].y;
		}
		for(int i=0; i < 1000; i++)
			ret = Solve5PointEssential(pts1, pts2, num_pts, Evec, Pvec, inliers);

		if(ret) {
			int maxsol=1,maxinliers=0;
			for (int i=0; i<inliers.size(); i++){
				if(inliers[i]>maxinliers){
					maxinliers=inliers[i];
					maxsol=i;
				}
			}
			R = Pvec[maxsol](cv::Rect(0,0,3,3));

			//cout << "Success! Found " <<  Evec.size() << " possible solutions" << endl;
			//cout << "The best one has the highest inliers. An inlier is a point that is in front of both cameras." << endl;
			//cout << endl;

			//for(size_t i=0; i < Evec.size(); i++) {
			//	cout << "Solution " << (i+1) << "/" << E.size() << endl;
			//	cout << endl;
			//	cout << "E = " << Evec[i] << endl;
			//	cout << endl;

			//	if(determinant(Pvec[i](Range(0,3), Range(0,3))) < 0) {
			//		cout << "Detected a reflection in P. Multiplying all values by -1." << endl;
			//		cout << "P = " << (Pvec[i] * -1) << endl;
			//	}
			//	else {
			//		cout << "P = " << Pvec[i] << endl;
			//	}

			//	cout << endl;
			//	cout << "inliers = " << inliers[i] << "/" << num_pts << endl;
			//	cout << "=========================================" << endl;
			//	cout << endl;
			//}
		}
		else {
			cout << "Could not find a valid essential matrix" << endl;
		}
 */


//stats.matches = (int)matched1.size();

//Mat inlier_mask, homography;
//vector<KeyPoint> inliers1, inliers2;
//vector<DMatch> inlier_matches;
//if(matched1.size() >= 4) {
//    homography = findHomography(Points(matched1), Points(matched2),
//                                RANSAC, ransac_thresh, inlier_mask);
//}

//if(matched1.size() < 4 || homography.empty()) {
//    Mat res;
//    hconcat(first_frame, frame, res);
//    stats.inliers = 0;
//    stats.ratio = 0;
//    return res;
//}
//for(unsigned i = 0; i < matched1.size(); i++) {
//    if(inlier_mask.at<uchar>(i)) {
//        int new_i = static_cast<int>(inliers1.size());
//        inliers1.push_back(matched1[i]);
//        inliers2.push_back(matched2[i]);
//        inlier_matches.push_back(DMatch(new_i, new_i, 0));
//    }
//}
//stats.inliers = (int)inliers1.size();
//stats.ratio = stats.inliers * 1.0 / stats.matches;

//vector<Point2f> new_bb;
//perspectiveTransform(object_bb, new_bb, homography);
//Mat frame_with_bb = frame.clone();
//if(stats.inliers >= bb_min_inliers) {
//    drawBoundingBox(frame_with_bb, new_bb);
//}
//Mat res;
//drawMatches(first_frame, inliers1, frame_with_bb, inliers2,
//            inlier_matches, res,
//            Scalar(255, 0, 0), Scalar(255, 0, 0));
//return res;
