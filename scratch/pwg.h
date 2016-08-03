#include "libs.h"
using namespace std;

class pwg {
    
public:
    
    struct point_3d {
        double x,y,z;
        point_3d(double x,double y,double z):x(x),y(y),z(z){}};
        
        struct point_2d {
            double x,y;
            point_2d(double x,double y):x(x),y(y){}};
            
            struct system {
                SpMat G,R,Y,H;
                VecXd x,y;
                vector<int> sw;} filter;
                
                vector<point_3d> triangulate(double *R, double *t, double *p1, double *p2, mwSize ncols);
                void formulate_system(vector<double> xc, vector<point_3d> xf);
                VecXd observation_model(mwSize ncols);
                SpMat observation_jacobian(mwSize ncols);
                void constraints_addition(VecXd z, VecXd zhat, SpMat H);
                void constraints_removal(VecXd z, VecXd zhat, SpMat H);
                
private:
    
    pair<VecXd,SpMat> recover_moments(VecXd y, SpMat Y);
    VecXd recover_first_moment(VecXd y, SpMat Y);
    SpMat recover_second_moment(SpMat Y);
    
    struct parameters{
        double sigma_s,sigma_b,sigma_a,sigma_r;
        double t1,t2;
        int exerror,debug,robust;
        double innov_gate,resid_gate;
        
        parameters (){
            sigma_s = pow (0.0050, 2);
            sigma_b = pow (0.0010, 2);
            sigma_a = pow (0.0087, 2); // 0.0087
            sigma_r = pow (0.0012, 2);
            t1      = 0.0001;
            t2      = 2.2204e-16;
            exerror = 0;
            debug   = 1;
            robust  = 1;
            innov_gate = 3.8415; //T = chi2inv(.95,1)
            resid_gate = 3.8415;}}; //T = chi2inv(.95,1)};
