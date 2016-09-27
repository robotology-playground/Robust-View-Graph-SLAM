#include <cpps.h>
#include <yarps.h>
#include <icubs.h>



//#include "Image.h"
//#include "ImageViewer.h"
#include "Tracker.h"
#include "PwgOptimiser.h"
#include "GraphOptimiser.h"
#include "featureselector.h"
#include "mythread.h"
//#include "5point.cpp"
//#include "Rpoly.cpp"

# define	M_PI	3.14159265358979323846  /* pi */

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
Network yarpnet; // 64 bytes still reachable
int i=0;

yarp::os::Mutex some_mutex;
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

bool mainThread(Mat& im,Tracker& t){
    t.setFirstFrame(im);
    return true;
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
     *          Then, generatgit commit -m "finished detector configuration part with ini files"
e_point_tracks should be integrated into
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

/* main code */
int main (int argc, char** argv) {
	/* we need at least one input, ncams */
	// argv[0] is the program name
    // argv[1:n] are the program input arguments
    int numImgs=20;
    double startTime=Time::now();
	if (argc<2) {
        std::cerr << "Usage: ./vgSLAM (int)ncams int(versbose)" << std::endl;
		return 1;
	}
	int ncams = std::atoi ( argv[1] );
    int VERBOSE = 0;
    if(ncams !=2 && ncams !=4){
        cout<<"ncams must be two or four"<<endl;
        return 1;
    }
//    ResourceFinder rf,rfdet,rfdesc,rfmatch;
//    if(!rf.setDefaultConfigFile("../../conf/vgSLAM.ini"))
//        return -1;
//    rf.configure(argc, argv);
//    rfdesc.configure(argc, argv);
//    rfmatch.configure(argc, argv);

    if (argc > 2){
        VERBOSE = std::atoi ( argv[2] );
    }

    FeatureSelector fs(argc,argv);
    Ptr<Feature2D> detector, descriptor;
    Ptr<DescriptorMatcher> matcher;

    fs.process(detector,descriptor,matcher);

	//if (VERBOSE==2){
	//	cvNamedWindow( "img_L", CV_WINDOW_AUTOSIZE );
	//	cvNamedWindow( "img_R", CV_WINDOW_AUTOSIZE );
	//}
    bool res=false;
	bool first=true;
    vector<Mat> imgvec(numImgs);
    while(i<numImgs){
        BufferedPort<ImageOf<PixelRgb> > image1_port, image2_port;
        int image1_start, image2_start;
        image1_port.open("/vgSLAM/cam/left");
        image2_port.open("/vgSLAM/cam/right");
        yarpnet.connect("/icub/cam/left", image1_port.getName());
        yarpnet.connect("/icub/cam/right",image2_port.getName());

        // initialise tracker
        std::vector<Tracker::point_2d> p ;

        //int nimages = 2*5063;
         //i<nimages

            ImageOf<PixelRgb> *image1_yarp = image1_port.read();
            ImageOf<PixelRgb> *image2_yarp = image2_port.read();
            Stamp s1,s2;
            if(image1_port.getEnvelope(s1) && image2_port.getEnvelope(s2)){
                if(i==0){
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
                    imgvec[i]=image1_cv;
                    imgvec[i+1]=image2_cv;
                    i=i+2;
                }
            }
    }

    std::cout<<"images acquired, now I analyse them in parallel"<<endl;
    Tracker tracker(detector,descriptor,matcher);
    int PrMatSize=0;
    if(ncams==4)
        PrMatSize=3*(numImgs-3);
    else
        PrMatSize=(numImgs-1);
    vector<Mat> ProjectionMatrices(PrMatSize);//Fix the size
    if(ncams==4){
        for(int j=0,i=0;j<numImgs-3;j++, i=i+3){
            yInfo()<<"Matching "<<j<<" to "<<j+1<<","<<j+2<<","<<j+3;
            cout<<"Waiting the first thread"<<endl;
            while(!mainThread(imgvec[j],tracker))
            {
                //do nothing.
            }
            MyThread t1(tracker,imgvec[j+1],ProjectionMatrices[i],1);//Fix the index
            MyThread t2(tracker,imgvec[j+2],ProjectionMatrices[i+1],2);//Fix the index
            MyThread t3(tracker,imgvec[j+3],ProjectionMatrices[i+2],3);//Fix the index
            t1.start();t2.start();t3.start();
            t1.stop();t2.stop();t3.stop();
            //cout<<ProjectionMatrices[j]<<endl;

        }
    }
    else {
        for(int j=0;j<numImgs-1;j++){
            yInfo()<<"Matching "<<j<<" to "<<j+1;
            while(!mainThread(imgvec[j],tracker))
            {
                //do nothing.
            }
            MyThread t1(tracker,imgvec[j+1],ProjectionMatrices[j],1);//FIX the management in the case of ncams==2;
            t1.start();
            t1.stop();
        }
    }

//    cout<<"Size"<<ProjectionMatrices.size()<<endl;
//    for(int i=0;i<ProjectionMatrices.size();i++)// OK TESTED
//        cout<<ProjectionMatrices[i]<<endl;
    double endTime=Time::now();
    std::cout <<endTime-startTime<< " seconds" << std::endl;
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
