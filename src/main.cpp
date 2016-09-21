#include <cpps.h>
#include <yarps.h>
#include <yarp/os/LogStream.h>
#include <icubs.h>
//#include <thread>
//#include <atomic>
//#include <mutex>
//#include <ctime>


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
//std::mutex some_mutex;
yarp::os::Mutex some_mutex;
/*
 *you may use Valgrind to test for any memory leaks
 *Valgrind is available for download at http://valgrind.org/
 *after compiling and installation, run
 *          valgrind ./PwgOptimiser 2 2
 *///for(int j=0;j<ncams;j=j+2){
//    cout<<"Serial version"<<endl;
//    Process2Images(detector,descriptor,matcher,imgvec[j],imgvec[j+1]);
//}
//vector<thread> threads;
//        for(int j=0; j < ncams; j=j+2 ){
//            thread t(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j]),ref(imgvec[j+1]),j/2);
////            cout << "main() : creating thread, " << j <<" with id: "<<t.get_id()<< endl;// qui l'id e' lo stesso. Quindi vuol dire che il thread e' uno e quindi no parallelismo.
//            if(t.joinable())
//                t.join();
//     threads.push_back(move(t));
//  }

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

bool Process2Images(Ptr<Feature2D>& detector, Ptr<Feature2D>& descriptor, Ptr<DescriptorMatcher>& matcher, Mat& image1_cv, Mat& image2_cv, int num)
{

                /* the tracker */
                cout<<"Starting the thread "<<num<<endl;
                Tracker tracker(detector, descriptor, matcher); // a tracker for each key-frame
                tracker.setFirstFrame(image1_cv);
                tracker.process(image2_cv);
                cout<<"Ending the thread "<<num<<endl;

}
//OLD deprecated
bool acquireAndProcess2Images(Ptr<Feature2D> detector, Ptr<Feature2D> descriptor, Ptr<DescriptorMatcher> matcher, int j){
    // YARP
    BufferedPort<ImageOf<PixelRgb> > image1_port, image2_port;
    int image1_start, image2_start;
    image1_port.open("/vgSLAM/cam/left"+to_string(j));
    image2_port.open("/vgSLAM/cam/right"+to_string(j));
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
                return false;
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
                tracker.setFirstFrame(image1_cv);//qui
                //	first = false;
                //}
                //else
                //	tracker.process(image1_cv);
                tracker.process(image2_cv);

                //get_aligned_point_matches (p, ncams, image1_cv) ;// thread safe
                //get_aligned_point_matches (p, ncams, image2_cv) ;// thread safe

                //ncams = p.size() ;
                //int npts = p[0].x.size() ;

                //if (VERBOSE==2){
                    //cvShowImage( "img_L", image1_cv ); // Large memory leak
                    //cvShowImage( "img_R", image2_cv ); // Large memory leak
                    //cvWaitKey(1);                     // Large memory leak
                    //cvDestroyWindow( "img_L" );
                    //cvDestroyWindow( "img_R" );
                //}

                //cvReleaseImage( &cvImageL );
                //cvReleaseImage( &cvImageR );
            }
            some_mutex.lock();
            i = i + 2;
            some_mutex.unlock();
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
    double startTime=Time::now();
	if (argc<2) {
        std::cerr << "Usage: ./vgSLAM (int)ncams int(versbose)" << std::endl;
		return 1;
	}
	int ncams = std::atoi ( argv[1] );
    int VERBOSE = 0;
    if(ncams%2==1){
        cout<<"ncams can be only multiple of two"<<endl;
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
    vector<Mat> imgvec(ncams);
    while(i<ncams){
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

    std::cout<<"images acquired, now I analyse them 2by2 in parallel"<<endl;
    // if ncams is multiple of four, use all the cpu power
    vector<MyThread*> threads;
    int cores=4;
    if(ncams%(cores*2)==0){
        for(int j=0;j<ncams;j=(j+cores*2)){
            cout<<"eight by eight images"<<endl;
            for (int i = 0; i < cores*2; i=i+2){
                MyThread* t=new MyThread(detector,descriptor,matcher,imgvec[j+i],imgvec[j+i+1],i/2);
                t->start();
                threads.push_back(t);
            }
            for (int i = 0; i < cores; ++i)
            {
                threads[i]->stop();
                delete threads[i];
                threads.pop_back();

            }

//            thread t1(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j]),ref(imgvec[j+1]),1);
//            thread t2(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j+2]),ref(imgvec[j+3]),2);
//            thread t3(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j+4]),ref(imgvec[j+5]),3);
//            thread t4(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j+6]),ref(imgvec[j+7]),4);
//            t1.join();t2.join();
//            t3.join();t4.join();
        }
    }
    else if(ncams%cores==0){
        cout<<"four by four images"<<endl;
        for(int j=0;j<ncams;j=j+cores){
            for (int i = 0; i < cores; i=i+2){
                MyThread* t=new MyThread(detector,descriptor,matcher,imgvec[j+i],imgvec[j+i+1],i/2);
                t->start();
                threads.push_back(t);
            }
            for (int i = 0; i < cores/2 ;++i)
            {
                threads[i]->stop();
                delete threads[i];
                threads.pop_back();
            }
//            thread t1(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j]),ref(imgvec[j+1]),1);
//            thread t2(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[j+2]),ref(imgvec[j+3]),2);
//            t1.join();t2.join();
        }
    }
    else
    {
        cout<<"two by two images"<<endl;
        for(int j=0;j<ncams;j=j+2){
            MyThread t1(detector,descriptor,matcher,imgvec[j],imgvec[j+1],0);
            t1.start();
            t1.stop();
        }
    }

//    thread t1(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[0]),ref(imgvec[1]),1);
//    //Process2Images(detector,descriptor,matcher,imgvec[2],imgvec[3],0); // neither main thread and child thread
//    thread t2(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[2]),ref(imgvec[3]),2);
//    //cout<<"hardware concurrency"<<t1.hardware_concurrency()<<end;
//    //cout<<"t1id: "<<t1.get_id()<<" t2id: "<<t2.get_id()<<endl; //e' diverso
//    thread t3(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[4]),ref(imgvec[5]),3);
//    thread t4(Process2Images,ref(detector),ref(descriptor),ref(matcher),ref(imgvec[6]),ref(imgvec[7]),4);
//    t1.join();t2.join();
//    t3.join();t4.join();
//for(int j=0;j<ncams;j=j+2){
//    cout<<"Serial version"<<endl;
//    Process2Images(detector,descriptor,matcher,imgvec[j],imgvec[j+1],j/2);
//}

	// run the process
	//process( ncams ) ;
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
            cout << "Could not find a valid essential matrix" << endl;            cout<<"QUI! threads size: "<<threads.size()<<endl;
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
