#ifndef PWGOPTIMISER_H
#define PWGOPTIMISER_H

#include "Eigen/Sparse"
#include "Eigen/Dense"
typedef Eigen::SparseMatrix<double> SpMat ;
typedef Eigen::Triplet<double> T ;

#include "RecoverMoments.h"

#include <unistd.h> // sleep command

class PwgOptimiser {
    
private:
    
    int ncams, npts ;
    Eigen::SparseMatrix<double> Y ;
    Eigen::VectorXd y ;

    struct constraint {
        const int cam ; /* camera identification number */
        const int kpt ; /* point identification number */
        const std::vector<double> p1 ; /* match in the first image */
        const std::vector<double> z ;  /* measurements or match in the second image */
        const std::vector<double> R ;  /* measurements noise covariance */
        Eigen::MatrixXd Y ;  /* measurements information matrix */
        Eigen::VectorXd y ;  /* measurements information vector */
        int sw; /* ON/OFF switch */
        constraint() ; // Do not use this default constructor (structure will be empty)
        constraint(
                int a, int b,
                std::vector<double> c,
                std::vector<double> d,
                std::vector<double> e,
                Eigen::MatrixXd F,
                Eigen::VectorXd G,
                int h):
                    cam(a), kpt(b), p1(c), z(d), R(e), Y(F), y(G) , sw(h) {} 
                    // Member initialization in a constructor
    };
    
    std::vector<constraint> constraints;
    //
    //std::vector<int> ON ;
    //std::vector<int> OFF ;
    std::vector<int> GATE_SWITCH ;
    std::vector<int> UPDATE_SWITCH ;
    //
    void initialise_info_matrix( 
            const double *xs );
    void constraints_addition_inverse_depth_Mviews( 
            const double *xs ) ;
    bool constraints_subtraction_inverse_depth_Mviews( 
            const double *xs ) ;
    void recover_moments ( 
            double *x, double *c ) ;
    void update_info_matrix_Mviews(
            Eigen::SparseMatrix<double>& Yon, Eigen::VectorXd& yon ) ;
    void compute_gate_inverse_depth_Mviews( //deals with Y elements
            double *gate, const double *x, const double *s, const double *xs ) ;
    void compute_gate_inverse_depth( //deals with Y as a sparse matrix
            double *gate, double *x, long unsigned int *ir,
            long unsigned int *jc, double *s, double *xs,
            long unsigned int ncol ) ;
    //
    void diagonal_triplet_form( 
            std::vector<T> &tripletList, int k) ;
    void off_diagonal_triplet_form( 
            std::vector<T> &tripletList, int k ) ;
    void update_information_vector( 
            Eigen::VectorXd& y, int k ) ;
    //
    void observation_model_inverse_depth_Mviews(
            double *z, const double *xs, double *p, int i, int c ) ;
    void observation_model_jacobian_inverse_depth_Mviews(
            Eigen::SparseMatrix<double> &H, 
            const double *xs, double *p,
            int i, int c ) ;
    //
    void observation_model_inverse_depth(
            double *z, double *x, double *p, long unsigned int N ) ;
    void observation_model_jacobian_inverse_depth(
            Eigen::SparseMatrix<double> &H,
            double *x, double *p, long unsigned int N ) ;
    //
    void compute_observation_model_derivatives(
            std::vector<double> &hi, std::vector<double> &xi,
            std::vector<double> &pi ) ;
    void compute_observation_model_derivatives_Mviews(
            std::vector<double> &hi, std::vector<double> &x,
            double *p ) ;
    void compute_rotation_matrix_derivatives(
            std::vector<double> &R, std::vector<double> &dRa,
            std::vector<double> &dRb, std::vector<double> &dRc,
            std::vector<double> &x ) ;
    //
    void transform_to_relative(
            std::vector<double> &x, std::vector<double> &p ) ;
    void compute_rotation_angles(
            std::vector<double> &R, std::vector<double> &x ) ;
    void compute_rotation_matrix(
            std::vector<double> &R, std::vector<double> &x ) ;
    void multiply(
            std::vector<double> &C, std::vector<double> &A,
            std::vector<double> &B, int nrows, int ncomm, int ncols ) ;
    void transpose(
            std::vector<double> &At, std::vector<double> &A,
            int nrows, int ncols ) ;
    void display_matrix(
            std::vector<double> &At, int nrows, int ncols ) ;
    void pi_to_pi(
            double x ) ;
    
public:
    
    // constructor
    PwgOptimiser( ) ;
    PwgOptimiser( int M, int N ) ;
    // destructor
    ~PwgOptimiser( ) ;
    //
    Eigen::SparseMatrix<double> Phat ;
    Eigen::VectorXd xhat ;
    //
    // this bypasses initialise_info_matrix by setting the information 
    // matrix and the information vector directly
    void set_information_matrix_and_vector(
            Eigen::VectorXd& yin, Eigen::SparseMatrix<double>& Yin ) ;
    void get_information_matrix_and_vector(
            Eigen::VectorXd& yout, Eigen::SparseMatrix<double>& Yout ) ;
    void get_switch_vector( int *sw ) ;
    //
    void initialise_a_constraint(const int& cam, const int& kpt, 
            const std::vector<double>& p1, const std::vector<double>& z, 
            const std::vector<double>& R, Eigen::MatrixXd Y, Eigen::VectorXd y, 
            const int& sw) ;
    void generate_constraints_info_Mviews( 
            const double *xs ) ;
    void optimise_constraints_image_inverse_depth_Mviews( 
            const double *xs ) ;
    struct pulled_constraint { /* A structure to pull a constraint for read only */
        const int cam;/* camera identification number */
        const int kpt;/* point identification number */
        const std::vector<double> p1 ;/* match in the first image */
        const std::vector<double> z ;/* measurements or match in the second image */
        const std::vector<double> R ;/* measurements noise covariance */
        const Eigen::MatrixXd Y ;/* measurements information matrix */
        const Eigen::VectorXd y ;/* measurements information vector */
        pulled_constraint() ;// Do not use this default constructor (structure will be empty)
        pulled_constraint(
                const constraint C):
                    cam(C.cam), kpt(C.kpt), p1(C.p1), z(C.z), R(C.R), Y(C.Y), y(C.y) {}
    } ;
    void pull_constraints_Mviews( 
            std::vector<pulled_constraint> &C ) ;
    //
    /* Allows mex functions to access private members */
    friend void mexFunction( 
            int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] );
    
//     friend void mex_observation_model_inverse_depth_Mviews(
//             double *z, const double *xs, double *p, int i, int c ) ;
//     friend void mex_observation_model_jacobian_inverse_depth_Mviews(
//             Eigen::SparseMatrix<double> &H, 
//             const double *xs, double *p,
//             int i, int c ) ;
//     friend void mex_observation_model_inverse_depth(
//             double *z, double *x, double *p, long unsigned int N ) ;
//     friend void mex_observation_model_jacobian_inverse_depth(
//             Eigen::SparseMatrix<double> &H,
//             double *x, double *p, long unsigned int N ) ;
//     friend void mex_compute_gate_inverse_depth_Mviews( //deals with Y elements
//             double *gate, const double *x, const double *s, const double *xs ) ;
//     friend void mex_compute_gate_inverse_depth( //deals with Y as a sparse matrix
//             double *gate, double *x, long unsigned int *ir,
//             long unsigned int *jc, double *s, double *xs,
//             long unsigned int ncol ) ;
//     friend void mex_constraints_addition_inverse_depth_Mviews( 
//             const double *xs ) ;
//     friend void mex_constraints_subtraction_inverse_depth_Mviews( 
//             const double *xs ) ;
//     friend void mex_recover_moments ( 
//             double *x, double *c ) ;
//     friend void mex_update_info_matrix_Mviews(
//             Eigen::SparseMatrix<double>& Yon, Eigen::VectorXd& yon ) ;
    
};

#endif /* PWGOPTIMISER_H */

/* 
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */