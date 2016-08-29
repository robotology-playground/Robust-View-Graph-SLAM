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
<<<<<<< HEAD
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

/* main code */
int main (int argc, char** argv) {
	/* we need at least one input, ncams */
	// argv[0] is the program name
	// argv[1:n] are the program input arguments

	if (argc<2) {
		std::cerr << "Usage: ./vgSLAM (int)ncams int(versbose)" << std::endl;
		return 1;
	}
	int ncams = std::atoi ( argv[1] );
	int VERBOSE = 0;
	if (argc > 2)
		VERBOSE = std::atoi ( argv[2] );
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
	Ptr<AKAZE> detector = AKAZE::create(); // features detector
	detector->setThreshold(akaze_thresh);
	Ptr<xfeatures2d::SIFT> descriptor = xfeatures2d::SIFT::create(); // features descriptor
	//Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("FlannBased"); // features matcher

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
				Tracker akaze_tracker(detector, descriptor, matcher); // a tracker for each key-frame
				//if(first){
				akaze_tracker.setFirstFrame(image1_cv);
				//	first = false;
				//}
				//else
				//	akaze_tracker.process(image1_cv);
				akaze_tracker.process(image2_cv);

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
=======
    double f = (double) std::rand() / RAND_MAX ;
    return fMin + f * (fMax - fMin) ;
}

/* image points structure */
struct point_2d {
    const std::vector<double> x ; /* x image coordinate */
    const std::vector<double> y ; /* y image coordinate */
    const std::vector<int> status ; /* status */
    int ncams = 0; // number of cameras
    int npts = 0; // number of cameras
    point_2d ( ) ; // Do not use this default constructor (structure will be empty)
    // Member initialization in a constructor
    point_2d ( std::vector<double> a, std::vector<double> b,
            std::vector<int> c) : x(a), y(b), status(c) {
                ncams++; npts=x.size(); }
};

/* generate image points */
void generate_point_tracks (std::vector<point_2d>& impoints, int ncams ) {
    std::vector<double> x, y ;
    std::vector<int> status ;
    
    /* FAST */
    int npts = 100 ;
    for (int j=0; j<npts; j++){
        x.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised x points
        y.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised y points
        status.push_back( 1 ) ;
    }
    impoints.push_back(point_2d( x, y, status ));
    x.clear() ;
    y.clear() ;
    status.clear() ;
    
    /* Optical flow */
    npts = impoints[0].x.size() ;
    for (int i=1; i<ncams; i++){
        for (int j=0; j<npts; j++){
            x.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised x points
            y.push_back( fRand ( 0.01 , 1 ) ) ; // generate random normalised y points
            status.push_back( 1 ) ;
        }
        impoints.push_back(point_2d( x, y, status ));
        x.clear() ;
        y.clear() ;
        status.clear() ;
    }
    
}

/* generates linearisation point */
double* generate_linearisation_point(std::vector<point_2d> p, double *Pkin){
    int ncams, npts ;
    ncams = p.size() ;
    npts = p[0].x.size() ;
    double *xs;
    //Pkin = generate_iCub_kinematics ( ncams ) ; /* generate forward kinematics matrices */
    //rho = generate_inverse_depth ( ncams ) ; /* generate point inverse depth */
    int i = 0 ;
    xs = new double [6*ncams+npts]() ;
    for (i=6*1; i<6*ncams; i=i+12)                 /* stereo constraints */
        xs[i] = 0.068 ; // Pkin goes here
    for (i=6*2; i<6*ncams; i=i+12)              /* monocular constraints */
        xs[i] = fRand(0.03, 0.05) ; // Pkin goes here
    for (i=6*ncams; i<6*ncams+npts; i++)     /* inverse depth parameters */
        xs[i] = fRand(1, 2) ; // rho goes here
    return xs ;
}

/* generates bundle constraints */
void generate_constraints_image_inverse_depth_Mviews(PwgOptimiser *Object,
        std::vector<point_2d> p){
    int ncams, npts ;
    ncams = p.size() ;
    npts = p[0].x.size() ;
    int cam, kpt, sw=0;
    std::vector<double> p1(2), z(2), R(4,0.0) ;
    Eigen::MatrixXd yz = Eigen::MatrixXd::Zero(7,1) ;
    Eigen::MatrixXd Yz = Eigen::MatrixXd::Zero(7,7) ;
    for (int i=0; i<ncams; i++)
        for (int j=0; j<npts; j++) {
            if (p[i].status[j]==1) {
                cam = i+1; // this should take values in the range 1:ncams
                kpt = j+1; // this should take values in the range 1:npts
                p1[0] = p[0].x[j] ; p1[1] = p[0].y[j] ;
                z[0] = p[i].x[j] ; z[1] = p[i].y[j] ;
                R[0] = 1 ; R[3] = 1 ;
                Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz, sw) ;
            }
        }
}

/* run the vision process */
void process (int ncams) {
    
    /* generate point tracks
     * here goes:
     *      1- Features extration
     *      2- Tracking or matching
     *      3- Visibility analysis
     *      4- Temporary/Real-time VO/kinematics solution ?
     */
    std::vector<point_2d> p ;
    generate_point_tracks (p,ncams ) ; // thread safe
    ncams = p.size() ;
    int npts = p[0].x.size() ;
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
    double *xs, *Pkin ;
    xs = generate_linearisation_point ( p, Pkin ) ;
    //for (int i=0; i<6*ncams+npts; i++) std::cout << xs[i] << std::endl ;

    /* initialise a PwgOptimiser object */
    PwgOptimiser *Object ; // pointer initialisation
    Object = new PwgOptimiser ( ncams, npts ) ; // pointer initialisation
    
    /* generate constraints */
    generate_constraints_image_inverse_depth_Mviews( Object, p ) ;

    /* optimise constraints information */
    Object->optimise_constraints_image_inverse_depth_Mviews( xs ) ;

    /* free memory */
    delete[] xs ;
    delete[] Pkin ;
    delete Object ; // delete class pointer
}

/* main code */
int main (int argc, char** argv) {
    /* we need at least one input, ncams */
    // argv[0] is the program name
    // argv[1:n] are the program input arguments
    
    if (argc<2) {
        std::cerr << "Usage: ./PwgOptimiser (int)ncams int(versbose) (string)dir" << std::endl;
        return 1;
    }
    
    int ncams = std::stoi ( argv[1] );
    int VERBOSE = 0;
    if (argc>2)
        VERBOSE = std::stoi ( argv[2] );
    
    // stream the images using yarpdataviewer
    if (argc<4) { // no input directory
        int i=0;
        bool res=false;
        bool first=true;
        double start;
        Network yarp; // 64 bytes still reachable
        
        BufferedPort<ImageOf<PixelRgb> > imagePortL,imagePortR;
        int startCountR=0, startCountL=0;
        imagePortL.open("/PwgOptimiser/cam_left");
        imagePortR.open("/PwgOptimiser/cam_right");
        yarp.connect("/icub/cam/left",imagePortL.getName());
        yarp.connect("/icub/cam/right",imagePortR.getName());
        
        int nimages = 10;
        
        if (VERBOSE==2){
            //cvNamedWindow( "img_L", CV_WINDOW_AUTOSIZE );
            //cvNamedWindow( "img_R", CV_WINDOW_AUTOSIZE );
        }
        
        while(i<nimages){
            
            ImageOf<PixelRgb> *yarpImageL = imagePortL.read();
            ImageOf<PixelRgb> *yarpImageR = imagePortR.read();
            
            Stamp sL,sR;
            
            if(imagePortL.getEnvelope(sL) && imagePortR.getEnvelope(sR)){
                
                if(first){
                    start=sL.getTime();
                    first=false;
                    //std::cout<<"Reinitialize"<<std::endl;
                }
            
                if(i==0){
                    startCountL = sL.getCount();
                    startCountR = sR.getCount();
                }
            
                if(abs(sL.getCount()-startCountL)>2 || abs(sR.getCount()-startCountR)>2 || fabs((sL.getTime())-(sR.getTime()))>0.03){//0.03 is the half delta t
                    startCountL = sL.getCount();
                    startCountR = sR.getCount();
                    continue;
                }
                
                if (yarpImageL!=NULL && yarpImageR!=NULL){
                    Image iL(*yarpImageL);
                    //IplImage* cvImageL = (IplImage*)iL.getIplImage();
                    Image iR(*yarpImageR);
                    //IplImage* cvImageR = (IplImage*)iR.getIplImage();
                    
                    //cvCvtColor(cvImageL, cvImageL, CV_RGB2BGR);
                    //cvCvtColor(cvImageR, cvImageR, CV_RGB2BGR);
                    
                    if (VERBOSE==2){
                        //cvShowImage( "img_L", cvImageL ); // Large memory leak
                        //cvShowImage( "img_R", cvImageR ); // Large memory leak
                        //cvWaitKey(1);                     // Large memory leak
                        //cvDestroyWindow( "img_L" );
                        //cvDestroyWindow( "img_R" );
                    }
                    
                    //cvReleaseImage( &cvImageL );
                    //cvReleaseImage( &cvImageR );
                }
                
                std::cout << "[" << i << "," << i+1 << "]" << std::endl;
                i = i + 2;
            }
        } 
    }
    
    else { // read the images from a folder
        // sequence directory
        std::string folder = argv[2];
        std::string folder_1 = folder+"/left/";
        std::string folder_2 = folder+"/right/";
        
        // count number of left images in the folder
        DIR *dir;
        struct dirent *ent;
        int count = 0;
        if ((dir = opendir (folder_1.c_str())) != NULL) {
            while ((ent = readdir (dir)) != NULL) {
                //cout << ent->d_name << endl;
                count++;
            }
            closedir (dir);
        } else {
            /* could not open directory */
            perror ("");
            return EXIT_FAILURE;
        }
        
        if (VERBOSE==2){
           //cvNamedWindow( "img_1", CV_WINDOW_AUTOSIZE );
           //cvNamedWindow( "img_2", CV_WINDOW_AUTOSIZE );
        }
        
        for (int32_t i=0; i<10; i++) {
            // input file names
            char base_name[256]; sprintf(base_name,"%08d.ppm",i);
            std::string left_img_file_name  = folder_1 + base_name;
            std::string right_img_file_name  = folder_2 + base_name;
            //IplImage* img_1 = cvLoadImage( left_img_file_name.c_str() ); // 88 bytes still reachable
            //IplImage* img_2 = cvLoadImage( right_img_file_name.c_str() ); // 88 bytes still reachable
            
            //cv::Mat image;
            //image = cv::imread(left_img_file_name.c_str(), CV_LOAD_IMAGE_COLOR);
            
            if (VERBOSE==2){
               //cvShowImage( "img_1", img_1 );
               //cvShowImage( "img_2", img_2 );
               //cvWaitKey(0);
            }
            
            //cvReleaseImage( &img_1 );
            //cvReleaseImage( &img_2 );
            
        }
        
        if (VERBOSE==2){
            //cvDestroyWindow( "img_1" );
            //cvDestroyWindow( "img_2" );
        }
        
    }
    
    // run the process
    process( ncams ) ;
    
    return 0 ;
>>>>>>> 3ff1d03e82a761413a402226243e64442986aaec
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

<<<<<<< HEAD
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
=======
//     //int id_1, id_2 ;
//     //Image im1, im2 ;
//
//     if (argc>2) {
//         for (int i = 2; i < argc; ++i) {
//             if (std::string(argv[i]) == "-left") {
//                 if (i + 1 < argc) { // Make sure we aren't at the end of argv!
//                     i++ ;
//                     //im1 = readPPM( argv[i] );
//                     //id_1 = i;
//                     IplImage* img_1 = cvLoadImage( argv[i] );
//                     cvNamedWindow( "img_1", CV_WINDOW_AUTOSIZE );
//                     cvShowImage( "img_1", img_1 );
//                     cvWaitKey(0);
//                     cvReleaseImage( &img_1 );
//                     cvDestroyWindow( "img_1" );
//
//                 }
//                 else {
//                     std::cerr<<"-left option requires one argument."<<std::endl;
//                     return 1;
//                 }
//             }
//             else if (std::string(argv[i]) == "-right") {
//                 if (i + 1 < argc) { // Make sure we aren't at the end of argv!
//                     i++;
//                     //im2 = readPPM( argv[i] );
//                     //id_2 = i;
//                     IplImage* img_2 = cvLoadImage( argv[i] );
//                     cvNamedWindow( "img_2", CV_WINDOW_AUTOSIZE );
//                     cvShowImage( "img_2", img_2 );
//                     cvWaitKey(0);
//                     cvReleaseImage( &img_2 );
//                     cvDestroyWindow( "img_2" );
//                 }
//                 else {
//                     std::cerr<<"-right option requires one argument."<<std::endl;
//                     return 1;
//                 }
//             }
//         }
//         //ImageViewer(argv[id_1], argv[id_2], im1.h, im1.w) ;
//     }
//
//     else {
//         std::cout << std::endl ;
//         std::cout << "Process iCub images " << std::endl;
//     }
>>>>>>> 3ff1d03e82a761413a402226243e64442986aaec
