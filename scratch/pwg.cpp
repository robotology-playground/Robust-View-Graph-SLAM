#include "pwg.h"
#include "common.cpp"
using namespace Eigen;

vector<pwg::point_3d> pwg::triangulate(double *R, double *t, double *p1, double *p2, size_t ncols){
    
    /* create the output matrix */
    vector<point_3d> xf;
    
    double a, b, c, d, e, f, denom, s, h;
    double D, r1, r2, rim, z;
    double tx=t[0], ty=t[1], tz=t[2];
    double r00=R[0], r01=R[3], r02=R[6];
    double r10=R[1], r11=R[4], r12=R[7];
    double r20=R[2], r21=R[5], r22=R[8];
    double v10, v11, v12, v20, v21, v22;
    double x0, y0, x1, y1;
    
    for (int i=0;i<ncols;i++){
        x0 = p1[2*(i+1)-2];
        y0 = p1[2*(i+1)-1];
        x1 = p2[2*(i+1)-2];
        y1 = p2[2*(i+1)-1];

        // calculate lines
        v10 = x0;
        v11 = y0;
        v12 = 1;
        v20 = r00*x1+r01*y1+r02;
        v21 = r10*x1+r11*y1+r12;
        v22 = r20*x1+r21*y1+r22;
        
        // distance formulas based on Schneider pp 409-412
		  // shortest distance between the two lines
        a = v10*v10+v11*v11+v12*v12;
        b = v10*v20+v11*v21+v12*v22;
        c = v20*v20+v21*v21+v22*v22;
        d = tx*v10+ty*v11+tz*v12;
        e = tx*v20+ty*v21+tz*v22;
        f = tx*tx+ty*ty+tz*tz;
        denom = a*c-b*b;
        if (denom < std::numeric_limits<double>::epsilon()) {
            denom = 1;} // accounts for parallel lines
        s = (c*d-b*e)/denom;
        h = (b*d-a*e)/denom;
        
        D = f+s*(a*s-b*h-2*d)+h*(c*h-b*s+2*e);
        
        // compute lines length
        r1 = sqrt(a)*s;
        r2 = sqrt(b)*h;
        
        // negative distance
        //vis = r1 > 0;
        
        // 3d points in first camera coordinates
        rim = sqrt(x0*x0+y0*y0+1);
        z = r1/rim;
        xf.push_back(point_3d(x0*z, y0*z, z));}
    
    // return 3d point vector
    return xf;}

void pwg::formulate_system(vector<double> xc, vector<point_3d> xf){
    
    parameters param;
    int index, nz, i, j;
    size_t n = xf.size();
    
    /* state vector and measurements switch */
    VecXd x(3*n+6);
    index = 0;
    i = 0; j = 0;
    while (index < n){
        x.coeffRef(i)=xf[index].x; i++;
        x.coeffRef(i)=xf[index].y; i++;
        x.coeffRef(i)=xf[index].z; i++;
        index++;}
    index = 0;
    while (index < 6){
        x.coeffRef(i)=xc[index];
        i++;
        index++;}
    pwg::filter.x=x;
    
    /* information matrix */
    vector<T> y_tripletList; // number of non-zeros
    index = 0;
    i = 1;
    j = 1;
    while (index < 3*n+6){
        if (index < 3*n){
            y_tripletList.push_back(T(i-1,j-1,1/param.sigma_s));}
        else if (index < 3*n+3){
            y_tripletList.push_back(T(i-1,j-1,1/param.sigma_b));}
        else if (index < 3*n+6){
            y_tripletList.push_back(T(i-1,j-1,1/param.sigma_a));}
        index++; i++; j++;}
    SpMat Y(3*n+6,3*n+6);
    Y.setFromTriplets(y_tripletList.begin(), y_tripletList.end());
    pwg::filter.Y=Y;
    
    /* information vector */
    VecXd y(3*n+6);
    y=Y*x;
    pwg::filter.y=y;
    
    /* measurements noise */
    vector<T> r_tripletList;
    vector<T> g_tripletList;
    index = 0;
    i = 1;
    j = 1;
    while (index < 2*n){
        r_tripletList.push_back(T(i-1,j-1,param.sigma_r));
        g_tripletList.push_back(T(i-1,j-1,1/param.sigma_r));
        index++; i++; j++;}
    SpMat R(2*n,2*n);
    R.setFromTriplets(r_tripletList.begin(), r_tripletList.end());
    pwg::filter.R=R;
    SpMat G(2*n,2*n);
    G.setFromTriplets(g_tripletList.begin(), g_tripletList.end());
    pwg::filter.G=G;

    /* observation switch */
    std::vector<int> sw(2*n); // all zeroes
    pwg::filter.sw = sw;}

VecXd pwg::observation_model(size_t ncols){
    /* state vector */
    VecXd x(3*ncols+6);
    x = pwg::filter.x;
    double tx=x.coeffRef(3*ncols+0);   // translations
    double ty=x.coeffRef(3*ncols+1);
    double tz=x.coeffRef(3*ncols+2);
    double rx=x.coeffRef(3*ncols+3); // rotations
    double ry=x.coeffRef(3*ncols+4);
    double rz=x.coeffRef(3*ncols+5);
    /* rotation matrix */
    double* R = (double*)_mm_malloc(3*3*sizeof(double),64);
    double* dR = (double*)_mm_malloc(9*3*sizeof(double),64);
    double* om = (double*)_mm_malloc(3*1*sizeof(double),64);
    om[0] = rx; om[1] = ry; om[2] = rz;
    common transform;
    transform.rodrigues(R, dR, om);
    double r00=R[0], r01=R[1], r02=R[2]; //R is transposed here
    double r10=R[3], r11=R[4], r12=R[5];
    double r20=R[6], r21=R[7], r22=R[8];
    double num1, num2, denum1;
    double tmp1,tmp2,tmp3;
    VecXd pix(2*ncols);
    for (int i=0;i<ncols;i++){
        tmp1 = x.coeffRef(3*i+0)-tx;
        tmp2 = x.coeffRef(3*i+1)-ty;
        tmp3 = x.coeffRef(3*i+2)-tz;
        num1   = r00*tmp1+r01*tmp2+r02*tmp3;
        num2   = r10*tmp1+r11*tmp2+r12*tmp3;
        denum1 = r20*tmp1+r21*tmp2+r22*tmp3;
        pix.coeffRef(2*i)  = num1/denum1;
        pix.coeffRef(2*i+1)= num2/denum1;}
    return pix;}

SpMat pwg::observation_jacobian(size_t ncols){
    /* state vector */
    VecXd x(3*ncols+6);
    x = pwg::filter.x;
    double tx=x.coeffRef(3*ncols+0);   // translations
    double ty=x.coeffRef(3*ncols+1);
    double tz=x.coeffRef(3*ncols+2);
    double rx=x.coeffRef(3*ncols+3); // rotations
    double ry=x.coeffRef(3*ncols+4);
    double rz=x.coeffRef(3*ncols+5);
    /* rotation matrix */
    double* R = (double*)_mm_malloc(3*3*sizeof(double),64);
    double* dR = (double*)_mm_malloc(9*3*sizeof(double),64);
    double* om = (double*)_mm_malloc(3*1*sizeof(double),64);
    om[0] = rx; om[1] = ry; om[2] = rz;
    common transform;
    transform.rodrigues(R, dR, om);
    double r00=R[0], r01=R[1], r02=R[2]; //R is transposed here
    double r10=R[3], r11=R[4], r12=R[5];
    double r20=R[6], r21=R[7], r22=R[8];
    double da00=dR[0], da01=dR[1], da02=dR[2]; //dRa is transposed here
    double da10=dR[3], da11=dR[4], da12=dR[5];
    double da20=dR[6], da21=dR[7], da22=dR[8];
    double db00=dR[9], db01=dR[10], db02=dR[11]; //dRb is transposed here
    double db10=dR[12], db11=dR[13], db12=dR[14];
    double db20=dR[15], db21=dR[16], db22=dR[17];
    double dc00=dR[18], dc01=dR[19], dc02=dR[20]; //dRc is transposed here
    double dc10=dR[21], dc11=dR[22], dc12=dR[23];
    double dc20=dR[24], dc21=dR[25], dc22=dR[26];
    double num1, num2, denum1, denum2, dnum1, dnum2, dnum3;
    double dudx, dudy, dudz, dvdx, dvdy, dvdz; // translation terms
    double duda, dudb, dudc, dvda, dvdb, dvdc; // rotation terms
    double dudX, dudY, dudZ, dvdX, dvdY, dvdZ; // map terms
    double tmp1,tmp2,tmp3;
    
    // initialise a triple list to fill in the sparse matrix J
    vector<T> tripletList;
    int skip1=0,skip2;
    for (int i=0;i<ncols;i++){
        tmp1 = x.coeffRef(3*i+0)-tx;
        tmp2 = x.coeffRef(3*i+1)-ty;
        tmp3 = x.coeffRef(3*i+2)-tz;
        num1   = r00*tmp1+r01*tmp2+r02*tmp3;
        num2   = r10*tmp1+r11*tmp2+r12*tmp3;
        denum1 = r20*tmp1+r21*tmp2+r22*tmp3;
        denum2 = denum1*denum1;
        // translations
        dudx = -r00/denum1 + num1*r20/denum2;
        dvdx = -r10/denum1 + num2*r20/denum2;
        dudy = -r01/denum1 + num1*r21/denum2;
        dvdy = -r11/denum1 + num2*r21/denum2;
        dudz = -r02/denum1 + num1*r22/denum2;
        dvdz = -r12/denum1 + num2*r22/denum2;
        // rotations
        dnum1 = da00*tmp1+da01*tmp2+da02*tmp3;
        dnum2 = da10*tmp1+da11*tmp2+da12*tmp3;
        dnum3 = da20*tmp1+da21*tmp2+da22*tmp3;
        duda = dnum1/denum1 - num1*dnum3/denum2;
        dvda = dnum2/denum1 - num2*dnum3/denum2;
        dnum1 = db00*tmp1+db01*tmp2+db02*tmp3;
        dnum2 = db10*tmp1+db11*tmp2+db12*tmp3;
        dnum3 = db20*tmp1+db21*tmp2+db22*tmp3;
        dudb = dnum1/denum1 - num1*dnum3/denum2;
        dvdb = dnum2/denum1 - num2*dnum3/denum2;
        dnum1 = dc00*tmp1+dc01*tmp2+dc02*tmp3;
        dnum2 = dc10*tmp1+dc11*tmp2+dc12*tmp3;
        dnum3 = dc20*tmp1+dc21*tmp2+dc22*tmp3;
        dudc = dnum1/denum1 - num1*dnum3/denum2;
        dvdc = dnum2/denum1 - num2*dnum3/denum2;
        // X-3d point terms
        dudX = r00/denum1 - num1*r20/denum2;
        dvdX = r10/denum1 - num2*r20/denum2;
        // Y-3d point terms
        dudY = r01/denum1 - num1*r21/denum2;
        dvdY = r11/denum1 - num2*r21/denum2;
        // Z-3d point terms
        dudZ = r02/denum1 - num1*r22/denum2;
        dvdZ = r12/denum1 - num2*r22/denum2;
//           (3*i) (3*i+1) (3*i+2) ... (3*ncols+0) ( +1) ( +2) ( +3) ( +4) ( +5)
//          ---------------------------------------------------------------------
// (2*i  )  | dudX   dudY   dudZ           dudx    dudy  dudz  duda  dudb  dudc |
// (2*i+1)  | dvdX   dvdY   dvdZ           dvdx    dvdy  dvdz  dvda  dvdb  dvdc |
//          ---------------------------------------------------------------------
        // Jacobian - map
        tripletList.push_back( T( 2*i+0, 3*i+0, dudX ) );
        tripletList.push_back( T( 2*i+1, 3*i+0, dvdX ) );
        tripletList.push_back( T( 2*i+0, 3*i+1, dudY ) );
        tripletList.push_back( T( 2*i+1, 3*i+1, dvdY ) );
        tripletList.push_back( T( 2*i+0, 3*i+2, dudZ ) );
        tripletList.push_back( T( 2*i+1, 3*i+2, dvdZ ) );
        // Jacobian - translations
        tripletList.push_back( T( 2*i+0, 3*ncols+0, dudx) );
        tripletList.push_back( T( 2*i+1, 3*ncols+0, dvdx) );
        tripletList.push_back( T( 2*i+0, 3*ncols+1, dudy) );
        tripletList.push_back( T( 2*i+1, 3*ncols+1, dvdy) );
        tripletList.push_back( T( 2*i+0, 3*ncols+2, dudz) );
        tripletList.push_back( T( 2*i+1, 3*ncols+2, dvdz) );
        // Jacobian - rotations
        tripletList.push_back( T(2*i+0, 3*ncols+3, duda) );
        tripletList.push_back( T(2*i+1, 3*ncols+3, dvda) );
        tripletList.push_back( T(2*i+0, 3*ncols+4, dudb) );
        tripletList.push_back( T(2*i+1, 3*ncols+4, dvdb) );
        tripletList.push_back( T(2*i+0, 3*ncols+5, dudc) );
        tripletList.push_back( T(2*i+1, 3*ncols+5, dvdc) );}
    //Assign them to the sparse Eigen matrix
    SpMat J(2*ncols,3*ncols+6);
    J.setFromTriplets(tripletList.begin(), tripletList.end());
    return J;}

void pwg::constraints_addition(VecXd z, VecXd zhat, SpMat H){
    parameters param;
    size_t n = z.size()/2;
    int k;
    SpMat R  = pwg::filter.R;
    SpMat G  = pwg::filter.G;
    SpMat Y  = pwg::filter.Y;
    VecXd xs = pwg::filter.x; // linearisation point
    VecXd y  = pwg::filter.y;
    std::vector<int> sw = pwg::filter.sw;
    
    int num_inliers = std::accumulate(sw.begin(),sw.end(),0);
    std::cout << "Information addition with initial constraints : " << num_inliers<< "\t";
    
    /* recover the second moment only */
    SpMat P = recover_second_moment(Y);
    
    /* Compute the innovations of the constraints in {z} that are currently OFF.
     * That is, for each sw that is OFF, compute the NIS for z */
    VecXd v(H.innerSize());
    v = z - zhat + H*xs;
    SpMat Ht_ = H.transpose();
    for (k=0;k<Ht_.outerSize();k++){
        if(sw[k]==1){ // remove sw=ON
            for (SparseMatrix<double,ColMajor>::InnerIterator it(Ht_,k);it;++it){
                Ht_.coeffRef(it.row(),it.col())=0;}}}
    SpMat H_ = Ht_.transpose();
    SpMat S = H_*P*Ht_+R;
    
    /* For each innovation in NIS, if NIS_j < T, turn ON the associated switch
     * sw_j and add the information due to zj to the state estimate, */
    double nis;
    std::vector<int> on(H.innerSize());
    SpMat Ht = H.transpose();
    for (k=0;k<Ht.outerSize();k++){
        if (sw[k]==0){
            nis=v.coeffRef(k)*v.coeffRef(k)/S.coeffRef(k,k);
            if(nis<param.innov_gate){
                on[k] = 1; sw[k] = 1;} 
            else { on[k] = 0;}} 
        else { on[k] = 0;}
        if(on[k] == 0){
            for (SparseMatrix<double,ColMajor>::InnerIterator it(Ht,k);it;++it){
                Ht.coeffRef(it.row(),it.col())=0;}}}
    H = Ht.transpose();
    int num_gate = std::accumulate(on.begin(),on.end(),0);
    std::cout << "added : " << num_gate << std::endl;
    
    y = y + Ht*G*v;
    Y = Y + Ht*G*H;
    pwg::filter.y  = y;
    pwg::filter.Y  = Y;
    pwg::filter.sw = sw;}

void pwg::constraints_removal(VecXd z, VecXd zhat, SpMat H){
    parameters param;
    size_t n = z.size()/2;
    int k;
    SpMat Y  = pwg::filter.Y;
    SpMat R  = pwg::filter.R;
    SpMat G  = pwg::filter.G;
    VecXd y  = pwg::filter.y;
    VecXd xs = pwg::filter.x; // linearisation point
    std::vector<int> sw = pwg::filter.sw;
    
    int num_inliers = std::accumulate(sw.begin(),sw.end(),0);
    std::cout << "Information removal with constraints : " << num_inliers << "\t";
    
    std::pair<VecXd,SpMat> pair = recover_moments(y, Y);
    VecXd x = pair.first;
    SpMat P = pair.second;
    
    /* Compute the residuals of the constraints in {z} that are currently ON.
     * That is, for each sw that is ON, compute the NRS for z */
    SpMat H_ = H;
    VecXd r(H.innerSize());
    r = z - zhat - H_*(x-xs);
    SpMat Ht_ = H_.transpose();
    for (k=0;k<Ht_.outerSize();k++){
        if(sw[k]==0){ // remove measurement if sw=0
            for (SparseMatrix<double,ColMajor>::InnerIterator it(Ht_,k);it;++it){
                Ht_.coeffRef(it.row(),it.col())=0;}}}
    H_ = Ht_.transpose();
    SpMat E = H_*P*Ht_+R;
    
    /* For each residual in NRS, if NRS_j > T, turn OFF the associated switch
     * sw_j and subtract the information due to zj from the state estimate,*/
    double nrs;
    std::vector<int> off(H.innerSize());
    SpMat Ht = H.transpose();
    for (k=0;k<Ht.outerSize();k++){
        if (sw[k]==1){
            nrs=r.coeffRef(k)*r.coeffRef(k)/E.coeffRef(k,k);
            if(nrs>param.resid_gate){
                off[k] = 1; sw[k]  = 0;
            } else { off[k] = 0; }
        } else { off[k] = 0; }
        if(off[k] == 0){
            for (SparseMatrix<double,ColMajor>::InnerIterator it(Ht,k);it;++it){
                Ht.coeffRef(it.row(),it.col())=0;}}}
    H = Ht.transpose();
    num_inliers = std::accumulate(sw.begin(),sw.end(),0);
    int num_gate = std::accumulate(off.begin(),off.end(),0);
    std::cout << "removed : " << num_gate << std::endl;
    
    VecXd v = z - zhat + H*xs;
    y = y - Ht*G*v;
    Y = Y - Ht*G*H;
    
    /* recovering the first moment */
    x = recover_first_moment(y, Y);
    pwg::filter.x  = x; // move the linearisation point to new estimate
    pwg::filter.y  = y;
    pwg::filter.Y  = Y;
    pwg::filter.sw = sw;}

std::pair <VecXd,SpMat> pwg::recover_moments(VecXd y, SpMat Y){
    size_t n = Y.outerSize();
    SpMat I(n,n);
    I.setIdentity();
    //ConjugateGradient<SpMat> solver;
    //SparseQR<SpMat, COLAMDOrdering<int> > solver(Y);
    SparseLU<SpMat, COLAMDOrdering<int> > solver(Y);
    SpMat P(n,n);
    P = solver.solve(I);
    VecXd x = P*y;
    std::pair<VecXd,SpMat> pair;
    pair.first = x;
    pair.second = P;
    return pair;}

VecXd pwg::recover_first_moment(VecXd y, SpMat Y){
    //ConjugateGradient<SpMat> solver;
    //SparseQR<SpMat, COLAMDOrdering<int> > solver(Y);
    SparseLU<SpMat, COLAMDOrdering<int> > solver(Y);
    VecXd x = solver.solve(y);
    return x;}

SpMat pwg::recover_second_moment(SpMat Y){
    size_t n = Y.outerSize();
    SpMat I(n,n);
    I.setIdentity();
    //ConjugateGradient<SpMat> solver;
    //SparseQR<SpMat, COLAMDOrdering<int> > solver(Y);
    SparseLU<SpMat, COLAMDOrdering<int> > solver(Y);
    SpMat P(n,n);
    P = solver.solve(I);
    return P;}
