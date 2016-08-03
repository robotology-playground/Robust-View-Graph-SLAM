#include <fstream>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <iostream>
#include <stdint.h>
#include <iomanip>

#include <opencv2/opencv.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#include <dirent.h> // search a directory for files

#include <yarp/os/all.h>
//#include <yarp/os/Time.h>
#include <yarp/sig/all.h>
#include <yarp/sig/Image.h>
#include <yarp/sig/Vector.h>
#include <yarp/dev/PolyDriver.h>
#include <yarp/dev/IEncoders.h>
#include <yarp/math/Math.h>

#include <iCub/iKin/iKinFwd.h>
#include <iCub/iKin/iKinIpOpt.h>

//#include <Eigen/Dense>
//#include <Eigen/StdVector>

//#include "Image.h"
//#include "ImageViewer.h"
#include "PwgOptimiser.h"

using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;
using namespace iCub::iKin;
//using namespace cv;
//using namespace std;

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
std::vector<point_2d> generate_point_tracks ( int ncams ) {
    std::vector<point_2d> impoints;
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
    
    return impoints ;
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
    Eigen::VectorXd Yz = Eigen::MatrixXd::Zero(7,7) ;
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
    p = generate_point_tracks ( ncams ) ; // thread safe
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
    //process( ncams ) ;
    
    return 0 ;
}





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