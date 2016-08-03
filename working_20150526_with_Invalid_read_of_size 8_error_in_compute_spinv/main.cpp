#include <iostream> // std::cout
#include <string>

#include "PwgOptimiser.h"

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
    int cam, kpt;
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
                Object->initialise_a_constraint(cam, kpt, p1, z, R, Yz, yz) ;
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
    p = generate_point_tracks ( ncams ) ;
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
    delete Object ; // delete class pointer
    
}

/* main code */
int main (int argc, char** argv) {
    
    /* we need at least one input, ncams */
    // argv[0] is the program name
    // argv[1:n] are the program input arguments
    if (argc<2) {
        std::cerr << "Usage: ./PwgOptimiser (int)ncams" << std::endl;
        return 1;
    }
    int ncams = std::stoi ( argv[1] ) ;
    //std::cout << "Using bundles of " << ncams << " cameras." << std::endl;
    /* image file options, npts */
    if (argc>2) {
        for (int i = 2; i < argc; ++i) {
            if (std::string(argv[i]) == "-left") {
                if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                    //left = argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
                    std::cout << "Read left image: " << argv[++i] << std::endl ;
                }
                else { // Uh-oh, there was no argument to the destination option.
                    std::cerr << "-left option requires one argument." << std::endl ;
                    return 1;
                }
            }
            else if (std::string(argv[i]) == "-right") {
                if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                    //right = argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
                    std::cout << "Read right image: " << argv[++i] << std::endl ;
                }
                else { // Uh-oh, there was no argument to the destination option.
                    std::cerr << "-right option requires one argument." << std::endl ;
                    return 1;
                }
            }
        }
    }
    else {
        std::cout << std::endl ;
        std::cout << "Process iCub images " << std::endl ;
    }
    
    process( ncams ) ;
    
    return 0 ;
}