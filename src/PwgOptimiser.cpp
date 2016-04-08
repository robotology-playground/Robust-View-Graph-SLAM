#include "PwgOptimiser.h"

/*=====================*/
PwgOptimiser::PwgOptimiser(void){
    //std::cout << "Object is being created" << std::endl;
}

/*=====================*/
PwgOptimiser::~PwgOptimiser(void){
    //std::cout << "Object is being deleted" << std::endl;
}

/*=====================*/
void PwgOptimiser::initialise_a_constraint(
        const int& cam,
        const int& kpt,
        const std::vector<double>& p1,
        const std::vector<double>& z,
        const std::vector<double>& R,
        Eigen::MatrixXd Y,
        Eigen::VectorXd y){
    constraints.push_back(PwgOptimiser::constraint(cam, kpt, p1, z, R, Y, y));
}

/*=====================*/
void PwgOptimiser::pull_constraints_Mviews(std::vector<pulled_constraint> &C){
    for (int ii=0; ii < constraints.size(); ii++) // constraints is the private structure
        C.push_back(pulled_constraint(constraints[ii]));
}

/*=====================*/
void PwgOptimiser::generate_constraints_info_Mviews(
        double *xs, long unsigned int ncams){
    
    Eigen::VectorXd dx(7), v(2), y(7,1);
    Eigen::MatrixXd R(2,2), H(2,7), Y(7,7);
    Eigen::SparseMatrix<double> Hs(2,7);
    double *p1 = new double [2];
    double *zs = new double [2];
    int count = 0;
    for (int ii=0; ii < constraints.size(); ii++){
        int c = (constraints[ii].cam-1)*6;
        int i = constraints[ii].kpt+6*ncams-1; // start from 0 as a C++ index
        int b = 0; // numerical checking flag
        constraints[ii].y = Eigen::MatrixXd::Zero(7,1);
        constraints[ii].Y = Eigen::MatrixXd::Zero(7,7);
        // Numerical checks (remove zero inverse depth)
        if (xs[i]==0){
            count++; 
            continue;
        }
        // reference 2d point
        for (int j=0; j<2; j++)//{
            p1[j] = constraints[ii].p1[j];
        // measurement covariance
        for (int j=0; j<2; j++){
           for (int k=0; k<2; k++){
               R.coeffRef(k,j) = constraints[ii].R[j*2+k];
           }
        };// std::cout << '\n' << R << std::endl;
        // jacobian
        observation_model_jacobian_inverse_depth_Mviews(Hs, xs, p1, i, c);
        H = Eigen::MatrixXd(Hs);
        // predicted measurements
        observation_model_inverse_depth_Mviews(zs, xs, p1, i, c);
        // state distance
        for (int j=0; j<6; j++)
            dx.coeffRef(j) = -xs[c+j];
        dx.coeffRef(6) = -xs[i];
        // innovations
        v = H*dx;
        for (int k=0; k<2; k++)
           v.coeffRef(k) = constraints[ii].z[k] - (zs[k] + v.coeffRef(k));
        //std::cout << '\n' << v << std::endl;
        // information vector
        y = (H.transpose()*R.inverse())*v;
        //std::cout << '\n' << y << std::endl;
        // information matrix
        Y = (H.transpose()*R.inverse())*H;
        //std::cout << '\n' << Y << std::endl;
        
        // Numerical checks (remove large information)
        for (int j=0;j<7;j++){
            if (Y.coeffRef(j,j)>1e+9)
                b = 1; 
        }
        if (b==1){
            count++; 
            continue;
        }
        // update constrain information
        constraints[ii].y = y;
        constraints[ii].Y = Y;
    }
    if (count>1)
        std::cout << count << " constraints were suppressed\n";
    else if (count==1)
        std::cout << count << " constraint was suppressed\n";
    delete[] p1, zs;
    
}

/*=====================*/
void PwgOptimiser::compute_gate_inverse_depth_Mviews(double *gate, 
        double *x, double *s, double *xs, long unsigned int ncols, 
        long unsigned int ncams){
    
    Eigen::VectorXd dx(7), v(2);
    Eigen::MatrixXd P(7,7), R(2,2), S(2,2), H(2,7);
    Eigen::SparseMatrix<double> Hs(2,7);
    double *p1 = new double [2];
    double *zs = new double [2];
    for (int ii=0; ii < constraints.size(); ii++){
        int c = (constraints[ii].cam-1)*6;
        int i = constraints[ii].kpt+6*ncams-1; // start from 0 as a C++ index
        p1[0] = constraints[ii].p1[0];
        p1[1] = constraints[ii].p1[1];
        // state distance
        for (int j=0; j<6; j++)
            dx.coeffRef(j) = x[c+j] - xs[c+j];
        dx.coeffRef(6) = x[i] - xs[i];
        pi_to_pi(dx.coeffRef(4));
        // state covariance
        for (int k=c; k<c+6; k++){
            for (int j=c; j<c+6; j++){
                P.coeffRef(j-c, k-c) = s[c*ncols + j + (k-c)*ncols];}
            P.coeffRef(k-c,6) = s[ncols*i + k];
            P.coeffRef(6,k-c) = s[ncols*i + k];}
        P.coeffRef(6,6) = s[ncols*i + i];
        //std::cout << '\n' << P << std::endl;
        // measurement covariance
        for (int j=0; j<2; j++){
           for (int k=0; k<2; k++){
               R.coeffRef(k,j) = constraints[ii].R[j*2+k];
           }
        };// std::cout << '\n' << R << std::endl;
        // jacobian
        observation_model_jacobian_inverse_depth_Mviews(Hs, xs, p1, i, c);
        H = Eigen::MatrixXd(Hs);
        //std::cout << '\n' << H << std::endl;
        // predicted measurements
        observation_model_inverse_depth_Mviews(zs, xs, p1, i, c);
        //std::cout << '\n';
        //for (int k=0; k<2; k++)
        //  std::cout << zs[k] << std::endl;
        // innovations
        v = H*dx;
        for (int k=0; k<2; k++)
           v.coeffRef(k) = constraints[ii].z[k] - (zs[k] + v.coeffRef(k));
        //std::cout << '\n' << v << std::endl;
        // innovations covariance
        S = H*P*H.transpose() + R;
        //std::cout << '\n' << S << std::endl;
        // gate function
        gate[ii] = v.transpose()*(S.inverse()*v);
    }
    delete[] p1, zs;
    
}

/*=====================*/
void PwgOptimiser::observation_model_inverse_depth_Mviews(double *z, double *xs,
        double *p, int i, int c){
    std::vector<double> pi(3,0.0), xi(6,0.0);
    double r, rim, d;
    r = 1/xs[i]; // range
    rim = sqrt(p[0]*p[0] + p[1]*p[1] + 1);
    d = r/rim;
    pi[0] = d*p[0]; // 3d-opints
    pi[1] = d*p[1];
    pi[2] = d;
    for (int i=0; i<6; i++)
            xi[i] = xs[c+i];
    transform_to_relative(xi, pi);
    z[0] = xi[0]/xi[2];
    z[1] = xi[1]/xi[2];
}

/*=====================*/
void PwgOptimiser::observation_model_jacobian_inverse_depth_Mviews(
        Eigen::SparseMatrix<double> &H,
        double *xs, double *p, int i, int c){
    std::vector<double> hi(14);
    std::vector<T> tripletList;
    std::vector<double> x(7,0.0); // 6 poses + 1 point
    for (int j=0; j<6; j++)
        x[j] = xs[c+j];
    x[6] = xs[i];
    compute_observation_model_derivatives_Mviews(hi, x, p);
    for (int j=0; j<7; j++){ // columns
        for (int i=0; i<2; i++){ // rows
            tripletList.push_back(T(i, j, hi[2*j+i]));
        }
    }
    H.setFromTriplets(tripletList.begin(), tripletList.end());
}

/*=====================*/
void PwgOptimiser::compute_observation_model_derivatives_Mviews(
        std::vector<double> &hi, std::vector<double> &x, double *p){

    // compute rotation matrix and its derivatives
    std::vector<double> x1(7,0.0); // 6 poses + 1 point
    for (int j=0; j<6; j++)
        x1[j] = x[j];
    std::vector<double> dRa(9,0.0),dRb(9,0.0),dRc(9,0.0),R1(9,0.0);
    compute_rotation_matrix_derivatives(R1,dRa,dRb,dRc,x1);
    
    // compute the derivatives
    double d, num1, num2, denum1, denum2, dnum1, dnum2, ddenum1, temp;
    std::vector<double> dp(3),xf(3);
    // 3d points in camera 1 frame
    d = 1/(x[6]*sqrt(p[0]*p[0]+p[1]*p[1]+1));
    dp[0] = -d*p[0]/x[6];
    dp[1] = -d*p[1]/x[6];
    dp[2] = -d*1/x[6];
    // Derivatives of camera 2 measurements
    double dx = d*p[0]-x[0];
    double dy = d*p[1]-x[1];
    double dz = d*1   -x[2];
    //display_matrix(dx,1,1);
    //display_matrix(dy,1,1);
    //display_matrix(dz,1,1);
    //display_matrix(R1,3,3);
    //display_matrix(R,3,3);
    num1   = R1[0]*dx + R1[1]*dy + R1[2]*dz;
    num2   = R1[3]*dx + R1[4]*dy + R1[5]*dz;
    denum1 = R1[6]*dx + R1[7]*dy + R1[8]*dz;
    denum2 = denum1*denum1;//denum2 = denum1.^2;
    dnum1   = R1[0]*dp[0] + R1[1]*dp[1] + R1[2]*dp[2];
    dnum2   = R1[3]*dp[0] + R1[4]*dp[1] + R1[5]*dp[2];
    ddenum1 = R1[6]*dp[0] + R1[7]*dp[1] + R1[8]*dp[2];
    // derivatives of pix with respect to translation t = [t1 t2 t3]'.
    hi[0] = -R1[0]/denum1 + num1*R1[6]/denum2; // dudt1
    hi[1] = -R1[3]/denum1 + num2*R1[6]/denum2; // dvdt1
    hi[2] = -R1[1]/denum1 + num1*R1[7]/denum2; // dudt2
    hi[3] = -R1[4]/denum1 + num2*R1[7]/denum2; // dvdt2
    hi[4] = -R1[2]/denum1 + num1*R1[8]/denum2; // dudt3
    hi[5] = -R1[5]/denum1 + num2*R1[8]/denum2; // dvdt3
    // derivatives of pix with respect to rotation a = [a1 a2 a3]'.
    xf[0] = dRa[0]*dx + dRa[3]*dy + dRa[6]*dz;
    xf[1] = dRa[1]*dx + dRa[4]*dy + dRa[7]*dz;
    xf[2] = dRa[2]*dx + dRa[5]*dy + dRa[8]*dz;
    hi[6] = xf[0]/denum1 - num1*xf[2]/denum2; // duda
    hi[7] = xf[1]/denum1 - num2*xf[2]/denum2; // dvda
    xf[0] = dRb[0]*dx + dRb[3]*dy + dRb[6]*dz;
    xf[1] = dRb[1]*dx + dRb[4]*dy + dRb[7]*dz;
    xf[2] = dRb[2]*dx + dRb[5]*dy + dRb[8]*dz;
    hi[8] = xf[0]/denum1 - num1*xf[2]/denum2; // dudb
    hi[9] = xf[1]/denum1 - num2*xf[2]/denum2; // dvdb
    xf[0] = dRc[0]*dx + dRc[3]*dy + dRc[6]*dz;
    xf[1] = dRc[1]*dx + dRc[4]*dy + dRc[7]*dz;
    xf[2] = dRc[2]*dx + dRc[5]*dy + dRc[8]*dz;
    hi[10] = xf[0]/denum1 - num1*xf[2]/denum2; // dudc
    hi[11] = xf[1]/denum1 - num2*xf[2]/denum2; // dvdc
    // derivatives of pix with respect to inverse depth rho = x(1)
    hi[12] = (dnum1*denum1 - num1*ddenum1)/denum2;
    hi[13] = (dnum2*denum1 - num2*ddenum1)/denum2;
}

/*=====================*/
void PwgOptimiser::compute_rotation_angles(std::vector<double> &R, std::vector<double> &x){
    x[3] = 0.5*(R[5]-R[7]);
    x[4] = 0.5*(R[6]-R[2]);
    x[5] = 0.5*(R[1]-R[3]);
    double theta = sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    if (theta>0.0000001){
        double trace = R[0] + R[4] + R[8];
        double w = atan2(theta,0.5*(trace-1));
        for (int i=3; i<6; i++)
            x[i] = x[i]/theta*w;
    }
}

/*=====================*/
void PwgOptimiser::compute_rotation_matrix(std::vector<double> &R, std::vector<double> &x){
    // compute rotation matrix and its derivatives
    double theta, alpha, beta, gamma;
    std::vector<double> omega(3),omegav(9),A(9);
    theta = sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    
    if (theta < 0.0000001){
        R[0] = 1; R[3] = 0; R[6] = 0;
        R[1] = 0; R[4] = 1; R[7] = 0;
        R[2] = 0; R[5] = 0; R[8] = 1;
    }
    else {
        omega[0] = x[3]/theta; omega[1] = x[4]/theta; omega[2] = x[5]/theta;
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1 - cos(theta);
        omegav[0] = 0;         omegav[3] = -omega[2];  omegav[6] = omega[1];
        omegav[1] = omega[2];  omegav[4] = 0;          omegav[7] = -omega[0];
        omegav[2] = -omega[1]; omegav[5] = omega[0];  omegav[8] = 0;
        A[0] = omega[0]*omega[0]; A[1] = omega[1]*omega[0];
        A[2] = omega[2]*omega[0];
        A[3] = omega[0]*omega[1]; A[4] = omega[1]*omega[1];
        A[5] = omega[2]*omega[1];
        A[6] = omega[0]*omega[2]; A[7] = omega[1]*omega[2];
        A[8] = omega[2]*omega[2];
        R[0] = alpha + omegav[0]*beta + A[0]*gamma;
        R[1] = omegav[1]*beta + A[1]*gamma;
        R[2] = omegav[2]*beta + A[2]*gamma;
        R[3] = omegav[3]*beta + A[3]*gamma;
        R[4] = alpha + omegav[4]*beta + A[4]*gamma;
        R[5] = omegav[5]*beta + A[5]*gamma;
        R[6] = omegav[6]*beta + A[6]*gamma;
        R[7] = omegav[7]*beta + A[7]*gamma;
        R[8] = alpha + omegav[8]*beta + A[8]*gamma;
    }
}

/*=====================*/
void PwgOptimiser::compute_rotation_matrix_derivatives(
        std::vector<double> &R,
        std::vector<double> &dRa,
        std::vector<double> &dRb,
        std::vector<double> &dRc,
        std::vector<double> &x){
    
    // compute rotation matrix and its derivatives
    double theta, alpha, beta, gamma;
    std::vector<double> omega(3),omegav(9),A(9);
    theta = sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    
    std::vector<double> dRdin(27,0.0);
    std::vector<double> dm3din(12,0.0);
    std::vector<double> dm2dm3(16,0.0);
    std::vector<double> dm1dm2(84,0.0);
    std::vector<double> dRdm1(189,0.0);
    std::vector<double> temp1(12,0.0);
    std::vector<double> temp2(63,0.0);
    
    if (theta<0.000000000000001){
        R[0] = 1; R[4] = 1; R[8] = 1;
        dRdin[5] = 1; dRdin[7] = -1; dRdin[11] = -1;
        dRdin[15] = 1; dRdin[19] = 1; dRdin[21] = -1;}
    else {
        omega[0] = x[3]/theta; omega[1] = x[4]/theta; omega[2] = x[5]/theta;
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1 - cos(theta);
        omegav[0] = 0;         omegav[3] = -omega[2];  omegav[6] = omega[1];
        omegav[1] = omega[2];  omegav[4] = 0;          omegav[7] = -omega[0];
        omegav[2] = -omega[1]; omegav[5] = omega[0];  omegav[8] = 0;
        A[0] = omega[0]*omega[0]; A[1] = omega[1]*omega[0];
        A[2] = omega[2]*omega[0];
        A[3] = omega[0]*omega[1]; A[4] = omega[1]*omega[1];
        A[5] = omega[2]*omega[1];
        A[6] = omega[0]*omega[2]; A[7] = omega[1]*omega[2];
        A[8] = omega[2]*omega[2];
        R[0] = alpha + omegav[0]*beta + A[0]*gamma;
        R[1] = omegav[1]*beta + A[1]*gamma;
        R[2] = omegav[2]*beta + A[2]*gamma;
        R[3] = omegav[3]*beta + A[3]*gamma;
        R[4] = alpha + omegav[4]*beta + A[4]*gamma;
        R[5] = omegav[5]*beta + A[5]*gamma;
        R[6] = omegav[6]*beta + A[6]*gamma;
        R[7] = omegav[7]*beta + A[7]*gamma;
        R[8] = alpha + omegav[8]*beta + A[8]*gamma;
        
        //std::vector<double> dm3din(12, 0.0);
        dm3din[0] = 1; dm3din[5] = 1; dm3din[10] = 1;
        dm3din[3] = omega[0]; dm3din[7] = omega[1]; dm3din[11] = omega[2];
        
        //std::vector<double> dm2dm3(16, 0.0);
        dm2dm3[0] = 1/theta; dm2dm3[12] = -omega[0]/theta;
        dm2dm3[5] = 1/theta; dm2dm3[13] = -omega[1]/theta;
        dm2dm3[10] = 1/theta; dm2dm3[14] = -omega[2]/theta;
        dm2dm3[15] = 1;
        
        //std::vector<double> dm1dm2(84,0.0); // initialised to all zeros (0.0)
        dm1dm2[63] = -sin(theta); dm1dm2[64] = cos(theta);
        dm1dm2[65] = sin(theta);
        dm1dm2[46] = 1; dm1dm2[26] = -1; dm1dm2[48] = -1;
        dm1dm2[8] = 1; dm1dm2[30] = 1; dm1dm2[10] = -1;
        dm1dm2[12] = 2*omega[0]; dm1dm2[13] = omega[1]; dm1dm2[34] = omega[0];
        dm1dm2[14] = omega[2]; dm1dm2[56] = omega[0];
        dm1dm2[15] = omega[1]; dm1dm2[36] = omega[0];
        dm1dm2[37] = 2*omega[1]; dm1dm2[38] = omega[2]; dm1dm2[59] = omega[1];
        dm1dm2[18] = omega[2]; dm1dm2[60] = omega[0];
        dm1dm2[40] = omega[2]; dm1dm2[61] = omega[1]; dm1dm2[62] = 2*omega[2];
        
        //std::vector<double> dRdm1(189,0.0); // initialised to all zeros (0.0)
        dRdm1[0] = 1; dRdm1[4] = 1; dRdm1[8] = 1;
        dRdm1[9]  = omegav[0]; dRdm1[10] = omegav[1]; dRdm1[11] = omegav[2];
        dRdm1[12] = omegav[3]; dRdm1[13] = omegav[4]; dRdm1[14] = omegav[5];
        dRdm1[15] = omegav[6]; dRdm1[16] = omegav[7]; dRdm1[17] = omegav[8];
        dRdm1[27] = beta; dRdm1[37] = beta;
        dRdm1[47] = beta; dRdm1[57] = beta;
        dRdm1[67] = beta; dRdm1[77] = beta;
        dRdm1[87] = beta; dRdm1[97] = beta;
        dRdm1[107] = beta;
        dRdm1[18] = A[0]; dRdm1[19] = A[1];
        dRdm1[20] = A[2]; dRdm1[21] = A[3];
        dRdm1[22] = A[4]; dRdm1[23] = A[5];
        dRdm1[24] = A[6]; dRdm1[25] = A[7];
        dRdm1[26] = A[8];
        dRdm1[108] = gamma; dRdm1[118] = gamma;
        dRdm1[128] = gamma; dRdm1[138] = gamma;
        dRdm1[148] = gamma; dRdm1[158] = gamma;
        dRdm1[168] = gamma; dRdm1[178] = gamma;
        dRdm1[188] = gamma;

        // dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;
        // 9x3     9x21     21x4     4x4      4x3
        // multiplies (4x4) x (4x3) = (21x3)
        multiply(temp1, dm2dm3, dm3din, 4, 4, 3);
        // multiplies (21x4) x (4x3) = (21x3)
        multiply(temp2, dm1dm2, temp1, 21, 4, 3);
        // multiplies (9x21) x (21x3) = (9x3)
        multiply(dRdin, dRdm1, temp2, 9, 21, 3);
    };
    // dRa
    dRa[0] = dRdin[0]; dRa[3] = dRdin[1]; dRa[6] = dRdin[2];
    dRa[1] = dRdin[3]; dRa[4] = dRdin[4]; dRa[7] = dRdin[5];
    dRa[2] = dRdin[6]; dRa[5] = dRdin[7]; dRa[8] = dRdin[8];
    
    // dRb
    dRb[0] = dRdin[9]; dRb[3] = dRdin[10]; dRb[6] = dRdin[11];
    dRb[1] = dRdin[12]; dRb[4] = dRdin[13]; dRb[7] = dRdin[14];
    dRb[2] = dRdin[15]; dRb[5] = dRdin[16]; dRb[8] = dRdin[17];
    
    // dRc
    dRc[0] = dRdin[18]; dRc[3] = dRdin[19]; dRc[6] = dRdin[20];
    dRc[1] = dRdin[21]; dRc[4] = dRdin[22]; dRc[7] = dRdin[23];
    dRc[2] = dRdin[24]; dRc[5] = dRdin[25]; dRc[8] = dRdin[26];
    
    // erase the vectors
    dRdin.erase(dRdin.begin(),dRdin.end());
    dm3din.erase(dm3din.begin(),dm3din.end());
    dm2dm3.erase(dm2dm3.begin(),dm2dm3.end());
    dm1dm2.erase(dm1dm2.begin(),dm1dm2.end());
    dRdm1.erase(dRdm1.begin(),dRdm1.end());
    temp1.erase(temp1.begin(),temp1.end());
    temp2.erase(temp2.begin(),temp2.end());
}

/*=====================*/
void PwgOptimiser::transform_to_relative(std::vector<double> &x, std::vector<double> &p){
    // compute rotation matrix
    double theta, alpha, beta, gamma;
    std::vector<double> omega(3), omegav(9,0.0), A(9), R(9);
    theta = sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    if (theta<0.000000000000001){
        R[0] = 1; R[4] = 1; R[8] = 1;}
    else {
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1 - cos(theta);
        omega[0] = x[3]/theta; omega[1] = x[4]/theta; omega[2] = x[5]/theta;
        omegav[3] = -omega[2]; omegav[6] =  omega[1]; omegav[7] = -omega[0];
        omegav[1] =  omega[2]; omegav[2] = -omega[1]; omegav[5] =  omega[0];
        A[0] = omega[0]*omega[0]; A[1] = omega[1]*omega[0]; A[2] = omega[2]*omega[0];
        A[3] = omega[0]*omega[1]; A[4] = omega[1]*omega[1]; A[5] = omega[2]*omega[1];
        A[6] = omega[0]*omega[2]; A[7] = omega[1]*omega[2]; A[8] = omega[2]*omega[2];
        R[0] = alpha + omegav[0]*beta + A[0]*gamma;
        R[1] = omegav[1]*beta + A[1]*gamma;
        R[2] = omegav[2]*beta + A[2]*gamma;
        R[3] = omegav[3]*beta + A[3]*gamma;
        R[4] = alpha + omegav[4]*beta + A[4]*gamma;
        R[5] = omegav[5]*beta + A[5]*gamma;
        R[6] = omegav[6]*beta + A[6]*gamma;
        R[7] = omegav[7]*beta + A[7]*gamma;
        R[8] = alpha + omegav[8]*beta + A[8]*gamma;}
    // compute global to relative transformation
    p[0] = p[0]-x[0];
    p[1] = p[1]-x[1];
    p[2] = p[2]-x[2];
    x[0] = R[0]*p[0]+R[1]*p[1]+R[2]*p[2];
    x[1] = R[3]*p[0]+R[4]*p[1]+R[5]*p[2];
    x[2] = R[6]*p[0]+R[7]*p[1]+R[8]*p[2];
}

/*=====================*/
void PwgOptimiser::multiply(std::vector<double> &C, std::vector<double> &A, std::vector<double> &B, int nrows, int ncomm, int ncols){
    // C = A*B
    // multiplies (nrows x ncomm) x (ncomm x ncols) = (nrows x ncols)
    for (int j = 0; j < ncols; j++){ // ncols of result
        for (int i = 0; i < nrows; i++){ // nrows of result
            for (int k = 0; k < ncomm; k++){ // common dimension
                C[nrows*j + i] = C[nrows*j + i]
                        + A[nrows*k + i] * B[ncomm*j + k];
            };
        };
    };
    for (int j = 0; j < nrows*ncols; j++){ // ncols of result
        if (sqrt(C[j]*C[j]) < 0.00000001)
                C[j] = 0;
    }
}

/*=====================*/
void PwgOptimiser::transpose(std::vector<double> &At, std::vector<double> &A,
        int nrows, int ncols){
    for (int j=0; j<ncols; j++){
        for (int i=0; i<nrows; i++){
            At[i*ncols+j] = A[ j*nrows+i];
        }
    }
}

/*=====================*/
void PwgOptimiser::display_matrix(std::vector<double> &A, int nrows, int ncols){
    std::cout << '\n';
    for (int j=0; j<nrows; j++){
        for (int i=0; i<ncols; i++){
            std::cout << A[i*nrows+j] << ' ';
        }
        std::cout << '\n';
    }
}

/*=====================*/
void PwgOptimiser::pi_to_pi(double &angle){
    double pi = 3.1415926534;
    double twopi = 2*pi;
    //angle = angle - twopi*fix(angle/twopi); // this is a stripped-down version of rem(angle, 2*pi)
    angle = angle - twopi*(angle - remainder(angle,twopi))/twopi;
    if (angle > pi){
        angle = angle - twopi;};
    if (angle < -pi){
        angle = angle + twopi;};
}

/*=====================*/
void PwgOptimiser::compute_gate_inverse_depth(double *gate, double *x,
        long unsigned int *ir, long unsigned int *jc, double *s, double *xs,
        long unsigned int ncol){
    
    long unsigned int ncon = constraints.size();
    
    Eigen::VectorXd dx(7), v(2), a(2), w(2);
    Eigen::MatrixXd P(7,7), R(2,2), S(2,2), H(2,7);
    Eigen::SparseMatrix<double> Hs(2,7);
    Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(7,7);
    R.coeffRef(0,0) = constraints[0].R[0];
    R.coeffRef(1,0) = 0.0;//constraints[0].R[1];
    R.coeffRef(0,1) = 0.0;//constraints[0].R[2];
    R.coeffRef(1,1) = constraints[0].R[0];
    
    //double *xi = (double *)calloc(7,sizeof(double));
    double *xi = new double [7];
    //double *pi = (double *)calloc(2,sizeof(double));
    double *pi = new double [2];
    //double *zi = (double *)calloc(2,sizeof(double));
    double *zi = new double [2];
    //double *zs = (double *)calloc(2,sizeof(double));
    double *zs = new double [2];
    xi[1] = xs[(ncol-6)+0];
    xi[2] = xs[(ncol-6)+1];
    xi[3] = xs[(ncol-6)+2];
    xi[4] = xs[(ncol-6)+3];
    xi[5] = xs[(ncol-6)+4];
    xi[6] = xs[(ncol-6)+5];
    dx.coeffRef(1) = x[(ncol-6)+0] - xs[(ncol-6)+0];
    dx.coeffRef(2) = x[(ncol-6)+1] - xs[(ncol-6)+1];
    dx.coeffRef(3) = x[(ncol-6)+2] - xs[(ncol-6)+2];
    dx.coeffRef(4) = x[(ncol-6)+3] - xs[(ncol-6)+3];
    dx.coeffRef(5) = x[(ncol-6)+4] - xs[(ncol-6)+4]; pi_to_pi(dx.coeffRef(4));
    dx.coeffRef(6) = x[(ncol-6)+5] - xs[(ncol-6)+5];
    
    for (int i=(ncol-6); i<(ncol-6)+6; i++) {            /* Loop through columns */
        for (int k=jc[i]; k<jc[i+1]; k++) {    /* Loop through non-zeros in ith column */
            if(ir[k]>(ncol-6)-1){ /* pick only pose elements */
                //std::cout << "row: " << ir[k]-(ncol-6)+1 << ", col: " << i-(ncol-6)+1 << ", data: " << s[k] << std::endl;
                P0.coeffRef(ir[k]-(ncol-6)+1, i-(ncol-6)+1) = s[k];}}}
    
    for (int ii=0; ii < ncon; ii++){
        
        int i = constraints[ii].kpt - 1; // fix to start from 0 as a C++ index
        xi[0] = xs[i];
        pi[0] = constraints[ii].p1[0];
        pi[1] = constraints[ii].p1[1];
        
        P = P0;
        for (int k=jc[i]; k<jc[i+1]; k++) {/* Loop through non-zeros in ith column */
            if(ir[k]==i){ /* pick only depth elements */
                //std::cout << "row: " << ir[k]-i << ", col: " << i-i << ", data: " << s[k] << std::endl;
                P.coeffRef(ir[k]-i, i-i) = s[k];}
            if(ir[k]>(ncol-6)-1){
                //std::cout << "row: " << ir[k]-(ncol-6)+1 << ", col: " << i-i << ", data: " << s[k] << std::endl;
                P.coeffRef(ir[k]-(ncol-6)+1, i-i) = s[k];}}
        for (int j=(ncol-6); j<(ncol-6)+6; j++) {/* Loop through columns */
            for (int k=jc[j]; k<jc[j+1]; k++) {/* Loop through non-zeros in ith column */
                if(ir[k]==i){ /* pick only depth elements */
                    //std::cout << "row: " << ir[k]-i << ", col: " << j-(ncol-6)+1 << ", data: " << s[k] << std::endl;
                    P.coeffRef(ir[k]-i, j-(ncol-6)+1) = s[k];}}}
        
        //Eigen::SparseMatrix<double> Hs(2,7);
        observation_model_jacobian_inverse_depth(Hs, xi, pi, 1);
        //Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2,7);
        H = Eigen::MatrixXd(Hs);
        observation_model_inverse_depth(zs, xi, pi, 1);
        
        dx.coeffRef(0) = x[i]-xs[i];
        a = H*dx;
        v.coeffRef(0) = constraints[ii].z[0] - zs[0] - a.coeffRef(0);
        v.coeffRef(1) = constraints[ii].z[1] - zs[1] - a.coeffRef(1);
        S = H*P*H.transpose()+R;
        
        //w = S.fullPivLu().solve(v); // works with singularity
        w = S.inverse()*v; // works with singularity
        //Eigen::LLT<Eigen::MatrixXd> llt;
        //llt.compute(S);
        //w = llt.solve(v);
        
        //double a = S.coeffRef(0);
        //double c = S.coeffRef(1);
        //double b = S.coeffRef(2);
        //double d = S.coeffRef(3);
        //S.coeffRef(0) = (1/(a*d-b*c))*d; // matrix inverse formulas
        //S.coeffRef(1) = -(1/(a*d-b*c))*c;
        //S.coeffRef(2) = -(1/(a*d-b*c))*b;
        //S.coeffRef(3) = (1/(a*d-b*c))*a;
        //w = S*v;
        
        gate[ii] = v.transpose()*w;
        //gate[ii] = v.coeffRef(0)*w.coeffRef(0)+v.coeffRef(1)*w.coeffRef(1);
        
    }
    
    delete[] xi, pi, zi, zs;
}

/*=====================*/
void PwgOptimiser::observation_model_jacobian_inverse_depth(
        Eigen::SparseMatrix<double> &H,
        double *x, double *p, long unsigned int N){
    std::vector<double> pi(3,0.0), xi(7,0.0), hi(14);
    std::vector<T> tripletList;
    xi[1] = x[N+0];
    xi[2] = x[N+1];
    xi[3] = x[N+2];
    xi[4] = x[N+3];
    xi[5] = x[N+4];
    xi[6] = x[N+5];
    for (int i = 0; i < N; i++){
        pi[0] = p[2*i+0];
        pi[1] = p[2*i+1];
        pi[2] = 1;
        xi[0] = x[i];
        compute_observation_model_derivatives(hi, xi, pi);
        tripletList.push_back(T(2*i+0, i, hi[0]));
        tripletList.push_back(T(2*i+1, i, hi[1]));
        tripletList.push_back(T(2*i+0, N+0, hi[2]));
        tripletList.push_back(T(2*i+1, N+0, hi[3]));
        tripletList.push_back(T(2*i+0, N+1, hi[4]));
        tripletList.push_back(T(2*i+1, N+1, hi[5]));
        tripletList.push_back(T(2*i+0, N+2, hi[6]));
        tripletList.push_back(T(2*i+1, N+2, hi[7]));
        tripletList.push_back(T(2*i+0, N+3, hi[8]));
        tripletList.push_back(T(2*i+1, N+3, hi[9]));
        tripletList.push_back(T(2*i+0, N+4, hi[10]));
        tripletList.push_back(T(2*i+1, N+4, hi[11]));
        tripletList.push_back(T(2*i+0, N+5, hi[12]));
        tripletList.push_back(T(2*i+1, N+5, hi[13]));}
    H.setFromTriplets(tripletList.begin(), tripletList.end());
}

/*=====================*/
void PwgOptimiser::observation_model_inverse_depth(double *z, double *x,
        double *p, long unsigned int N){
    std::vector<double> pi(3,0.0), xi(6,0.0);
    double r, rim, d, dd;
    for (int i=0; i<6; i++)
            xi[i] = x[N+i];
    for (int i = 0; i < N; i++){
        pi[0] = p[2*i+0];
        pi[1] = p[2*i+1];
        pi[2] = 1;
        r = 1/x[i]; // range
        rim = sqrt(pi[0]*pi[0] + pi[1]*pi[1] + 1);
        d = r/rim;
        pi[0] = d*pi[0]; // 3d-opints
        pi[1] = d*pi[1];
        pi[2] = d*pi[2];
        transform_to_relative(xi, pi);
        z[2*i+0] = xi[0]/xi[2];
        z[2*i+1] = xi[1]/xi[2];
    }
}

/*=====================*/
void PwgOptimiser::compute_observation_model_derivatives(
        std::vector<double> &hi,
        std::vector<double> &x,
        std::vector<double> &p){
    // compute rotation matrix and its derivatives
    double theta, alpha, beta, gamma;
    std::vector<double> dRa(9),dRb(9),dRc(9),omega(3),omegav(9),A(9),R(9);
    theta = sqrt(x[4]*x[4] + x[5]*x[5] + x[6]*x[6]);
    
    std::vector<double> dRdin(27,0.0);
    std::vector<double> dm3din(12,0.0);
    std::vector<double> dm2dm3(16,0.0);
    std::vector<double> dm1dm2(84,0.0);
    std::vector<double> dRdm1(189,0.0);
    std::vector<double> temp1(12,0.0);
    std::vector<double> temp2(63,0.0);
    
    if (theta<0.000000000000001){
        R[0] = 1; R[4] = 1; R[8] = 1;
        dRdin[5] = 1; dRdin[7] = -1; dRdin[11] = -1;
        dRdin[15] = 1; dRdin[19] = 1; dRdin[21] = -1;
    }
    else {
        omega[0] = x[4]/theta; omega[1] = x[5]/theta; omega[2] = x[6]/theta;
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1 - cos(theta);
        omegav[0] = 0;         omegav[3] = -omega[2];  omegav[6] = omega[1];
        omegav[1] = omega[2];  omegav[4] = 0;          omegav[7] = -omega[0];
        omegav[2] = -omega[1]; omegav[5] = omega[0];  omegav[8] = 0;
        A[0] = omega[0]*omega[0]; A[1] = omega[1]*omega[0];
        A[2] = omega[2]*omega[0];
        A[3] = omega[0]*omega[1]; A[4] = omega[1]*omega[1];
        A[5] = omega[2]*omega[1];
        A[6] = omega[0]*omega[2]; A[7] = omega[1]*omega[2];
        A[8] = omega[2]*omega[2];
        R[0] = alpha + omegav[0]*beta + A[0]*gamma;
        R[1] = omegav[1]*beta + A[1]*gamma;
        R[2] = omegav[2]*beta + A[2]*gamma;
        R[3] = omegav[3]*beta + A[3]*gamma;
        R[4] = alpha + omegav[4]*beta + A[4]*gamma;
        R[5] = omegav[5]*beta + A[5]*gamma;
        R[6] = omegav[6]*beta + A[6]*gamma;
        R[7] = omegav[7]*beta + A[7]*gamma;
        R[8] = alpha + omegav[8]*beta + A[8]*gamma;
        
        //std::vector<double> dm3din(12, 0.0);
        dm3din[0] = 1; dm3din[5] = 1; dm3din[10] = 1;
        dm3din[3] = omega[0]; dm3din[7] = omega[1]; dm3din[11] = omega[2];
        
        //std::vector<double> dm2dm3(16, 0.0);
        dm2dm3[0] = 1/theta; dm2dm3[12] = -omega[0]/theta;
        dm2dm3[5] = 1/theta; dm2dm3[13] = -omega[1]/theta;
        dm2dm3[10] = 1/theta; dm2dm3[14] = -omega[2]/theta;
        dm2dm3[15] = 1;
        
        //std::vector<double> dm1dm2(84,0.0); // initialised to all zeros (0.0)
        dm1dm2[63] = -sin(theta); dm1dm2[64] = cos(theta);
        dm1dm2[65] = sin(theta);
        dm1dm2[46] = 1; dm1dm2[26] = -1; dm1dm2[48] = -1;
        dm1dm2[8] = 1; dm1dm2[30] = 1; dm1dm2[10] = -1;
        dm1dm2[12] = 2*omega[0]; dm1dm2[13] = omega[1]; dm1dm2[34] = omega[0];
        dm1dm2[14] = omega[2]; dm1dm2[56] = omega[0];
        dm1dm2[15] = omega[1]; dm1dm2[36] = omega[0];
        dm1dm2[37] = 2*omega[1]; dm1dm2[38] = omega[2]; dm1dm2[59] = omega[1];
        dm1dm2[18] = omega[2]; dm1dm2[60] = omega[0];
        dm1dm2[40] = omega[2]; dm1dm2[61] = omega[1]; dm1dm2[62] = 2*omega[2];
        
        //std::vector<double> dRdm1(189,0.0); // initialised to all zeros (0.0)
        dRdm1[0] = 1; dRdm1[4] = 1; dRdm1[8] = 1;
        dRdm1[9]  = omegav[0]; dRdm1[10] = omegav[1]; dRdm1[11] = omegav[2];
        dRdm1[12] = omegav[3]; dRdm1[13] = omegav[4]; dRdm1[14] = omegav[5];
        dRdm1[15] = omegav[6]; dRdm1[16] = omegav[7]; dRdm1[17] = omegav[8];
        dRdm1[27] = beta; dRdm1[37] = beta;
        dRdm1[47] = beta; dRdm1[57] = beta;
        dRdm1[67] = beta; dRdm1[77] = beta;
        dRdm1[87] = beta; dRdm1[97] = beta;
        dRdm1[107] = beta;
        dRdm1[18] = A[0]; dRdm1[19] = A[1];
        dRdm1[20] = A[2]; dRdm1[21] = A[3];
        dRdm1[22] = A[4]; dRdm1[23] = A[5];
        dRdm1[24] = A[6]; dRdm1[25] = A[7];
        dRdm1[26] = A[8];
        dRdm1[108] = gamma; dRdm1[118] = gamma;
        dRdm1[128] = gamma; dRdm1[138] = gamma;
        dRdm1[148] = gamma; dRdm1[158] = gamma;
        dRdm1[168] = gamma; dRdm1[178] = gamma;
        dRdm1[188] = gamma;
        
        // dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;
        // 9x3     9x21     21x4     4x4      4x3
        //std::vector<double> temp1(12,0.0);
        temp1[0] = dm2dm3[0]*dm3din[0] + dm2dm3[4]*dm3din[1]
                + dm2dm3[8]*dm3din[2] + dm2dm3[12]*dm3din[3];
        temp1[1] = dm2dm3[1]*dm3din[0] + dm2dm3[5]*dm3din[1]
                + dm2dm3[9]*dm3din[2] + dm2dm3[13]*dm3din[3];
        temp1[2] = dm2dm3[2]*dm3din[0] + dm2dm3[6]*dm3din[1]
                + dm2dm3[10]*dm3din[2] + dm2dm3[14]*dm3din[3];
        temp1[3] = dm2dm3[3]*dm3din[0] + dm2dm3[7]*dm3din[1]
                + dm2dm3[11]*dm3din[2] + dm2dm3[15]*dm3din[3];
        temp1[4] = dm2dm3[0]*dm3din[4] + dm2dm3[4]*dm3din[5]
                + dm2dm3[8]*dm3din[6] + dm2dm3[12]*dm3din[7];
        temp1[5] = dm2dm3[1]*dm3din[4] + dm2dm3[5]*dm3din[5]
                + dm2dm3[9]*dm3din[6] + dm2dm3[13]*dm3din[7];
        temp1[6] = dm2dm3[2]*dm3din[4] + dm2dm3[6]*dm3din[5]
                + dm2dm3[10]*dm3din[6] + dm2dm3[14]*dm3din[7];
        temp1[7] = dm2dm3[3]*dm3din[4] + dm2dm3[7]*dm3din[5]
                + dm2dm3[11]*dm3din[6] + dm2dm3[15]*dm3din[7];
        temp1[8] = dm2dm3[0]*dm3din[8] + dm2dm3[4]*dm3din[9]
                + dm2dm3[8]*dm3din[10] + dm2dm3[12]*dm3din[11];
        temp1[9] = dm2dm3[1]*dm3din[8] + dm2dm3[5]*dm3din[9]
                + dm2dm3[9]*dm3din[10] + dm2dm3[13]*dm3din[11];
        temp1[10] = dm2dm3[2]*dm3din[8] + dm2dm3[6]*dm3din[9]
                + dm2dm3[10]*dm3din[10] + dm2dm3[14]*dm3din[11];
        temp1[11] = dm2dm3[3]*dm3din[8] + dm2dm3[7]*dm3din[9]
                + dm2dm3[11]*dm3din[10] + dm2dm3[15]*dm3din[11];
        
        // multiplies (21x4) x (4x3) = (21x3)
        int ncols = 3;
        int nrows = 21;
        int ncomm = 4;
        //std::vector<double> temp2(63,0.0);
        for (int j = 0; j < ncols; j++){ // ncols of result
            for (int i = 0; i < nrows; i++){ // nrows of result
                for (int k = 0; k < ncomm; k++){ // common dimension
                    temp2[nrows*j + i] = temp2[nrows*j + i]
                            + dm1dm2[nrows*k + i] * temp1[ncomm*j + k];
                };
            };
        };
        
        // multiplies (9x21) x (21x3) = (9x3)
        ncols = 3;
        nrows = 9;
        ncomm = 21;
        for (int j = 0; j < ncols; j++){ // ncols of result
            for (int i = 0; i < nrows; i++){ // nrows of result
                for (int k = 0; k < ncomm; k++){ // common dimension
                    dRdin[nrows*j + i] = dRdin[nrows*j + i]
                            + dRdm1[nrows*k + i] * temp2[ncomm*j + k];
                };
            };
        };
    };
    
    // compute the derivatives
    double r, rim, d, dd;
    double num1, num2, denum1, denum2, dnum1, dnum2, ddenum1, temp;
    std::vector<double> dp(3), xf(3);
    r = 1/x[0]; // range from inverse depth
    rim = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]); // 3d points in camera 1 frame
    d = r/rim;
    dd = -r*d;
    dp[0] = dd*p[0];
    dp[1] = dd*p[1];
    dp[2] = dd*p[2];
    p[0] = d*p[0];
    p[1] = d*p[1];
    p[2] = d*p[2];
    // Derivatives of camera 2 measurements
    p[0] = p[0] - x[1];
    p[1] = p[1] - x[2];
    p[2] = p[2] - x[3];
    num1 = R[0]*p[0] + R[1]*p[1] + R[2]*p[2];//num1 = R(1,:)*p;
    num2 = R[3]*p[0] + R[4]*p[1] + R[5]*p[2];//num2 = R(2,:)*p;
    denum1 = R[6]*p[0] + R[7]*p[1] + R[8]*p[2];//denum1 = R(3,:)*p;
    denum2 = denum1*denum1;//denum2 = denum1.^2;
    // R is transposed here:
    dnum1 = R[0]*dp[0] + R[1]*dp[1] + R[2]*dp[2];//dnum1 = R(1,:)*dp;
    dnum2 = R[3]*dp[0] + R[4]*dp[1] + R[5]*dp[2];//dnum2 = R(2,:)*dp;
    ddenum1 = R[6]*dp[0] + R[7]*dp[1] + R[8]*dp[2];//ddenum1 = R(3,:)*dp;
    // derivatives of pix with respect to inverse depth rho = x(1)
    hi[0] = (dnum1*denum1 - num1*ddenum1)/denum2;
    hi[1] = (dnum2*denum1 - num2*ddenum1)/denum2;
    // derivatives of pix with respect to translation t = [t1 t2 t3]'.
    temp = R[6]/denum2; // temp = R(3,1)./denum2;
    hi[2] = - R[0]/denum1 + num1*temp; // dudt1
    hi[3] = - R[3]/denum1 + num2*temp; // dvdt1
    temp = R[7]/denum2;
    hi[4] = - R[1]/denum1 + num1*temp; // dudt2
    hi[5] = - R[4]/denum1 + num2*temp; // dvdt2
    temp = R[8]/denum2;
    hi[6] = - R[2]/denum1 + num1*temp; // dudt3
    hi[7] = - R[5]/denum1 + num2*temp; // dvdt3
    // derivatives of pix with respect to rotation a = [a1 a2 a3]'.
    xf[0] = dRdin[0]*p[0] + dRdin[1]*p[1] + dRdin[2]*p[2];
    xf[1] = dRdin[3]*p[0] + dRdin[4]*p[1] + dRdin[5]*p[2];
    xf[2] = dRdin[6]*p[0] + dRdin[7]*p[1] + dRdin[8]*p[2];
    temp = xf[2]/denum2;
    hi[8]  = xf[0]/denum1 - num1*temp; // duda
    hi[9] = xf[1]/denum1 - num2*temp; // dvda
    xf[0] = dRdin[9]*p[0] + dRdin[10]*p[1] + dRdin[11]*p[2];
    xf[1] = dRdin[12]*p[0] + dRdin[13]*p[1] + dRdin[14]*p[2];
    xf[2] = dRdin[15]*p[0] + dRdin[16]*p[1] + dRdin[17]*p[2];
    temp = xf[2]/denum2;
    hi[10] = xf[0]/denum1 - num1*temp; // dudb
    hi[11] = xf[1]/denum1 - num2*temp; // dvdb
    xf[0] = dRdin[18]*p[0] + dRdin[19]*p[1] + dRdin[20]*p[2];
    xf[1] = dRdin[21]*p[0] + dRdin[22]*p[1] + dRdin[23]*p[2];
    xf[2] = dRdin[24]*p[0] + dRdin[25]*p[1] + dRdin[26]*p[2];
    temp = xf[2]/denum2;
    hi[12] = xf[0]/denum1 - num1*temp; // dudc
    hi[13] = xf[1]/denum1 - num2*temp; // dvdc
    
    // erase the vectors
    dRdin.erase(dRdin.begin(),dRdin.end());
    dm3din.erase(dm3din.begin(),dm3din.end());
    dm2dm3.erase(dm2dm3.begin(),dm2dm3.end());
    dm1dm2.erase(dm1dm2.begin(),dm1dm2.end());
    dRdm1.erase(dRdm1.begin(),dRdm1.end());
    temp1.erase(temp1.begin(),temp1.end());
    temp2.erase(temp2.begin(),temp2.end());
}