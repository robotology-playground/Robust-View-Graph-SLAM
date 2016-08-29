#include "GraphOptimiser.h"

/*=====================*/
GraphOptimiser::GraphOptimiser(void){
	//std::cout << "Object is being created" << std::endl;
}

/*=====================*/
GraphOptimiser::~GraphOptimiser(void){
	//std::cout << "Object is being deleted" << std::endl;
}

/*=====================*/
void GraphOptimiser::initialise_a_constraint(
		const std::vector<double>& edge,
		const std::vector<double>& z,
		const std::vector<double>& R){
	constraints.push_back(GraphOptimiser::constraint(edge, z, R));
}

/*=====================*/
void GraphOptimiser::compute_gate(double *gate, double *x, double *s,
		double *xs, long unsigned int ncols){

	long unsigned int ncons = constraints.size();
	Eigen::VectorXd dx(12), z(6), v(6);
	Eigen::MatrixXd P(12,12), R(6,6), S(6,6), H(6,12);
	Eigen::SparseMatrix<double> Hs(6,12);

	for (int ii=0; ii < ncons; ii++){

		// nodes
		int c1 = (constraints[ii].edge[0]-1)*6;
		int c2 = (constraints[ii].edge[1]-1)*6;

		// state distance
		for (int i=0; i<6; i++)
			dx.coeffRef(i) = x[c1+i] - xs[c1+i];
		pi_to_pi(dx.coeffRef(4));
		for (int i=0; i<6; i++)
			dx.coeffRef(i+6) = x[c2+i] - xs[c2+i];
		pi_to_pi(dx.coeffRef(10));
		//std::cout << '\n' << dx << std::endl;

		// state covariance
		for (int k=c1; k<c1+6; k++){ // defines the column
			for (int j=c1; j<c1+6; j++)// defines the row
				P.coeffRef(j-c1, k-c1) = s[k*ncols+j];
			for (int j=c2; j<c2+6; j++)
				P.coeffRef(j-c2+6, k-c1) = s[k*ncols+j];
		}
		for (int k=c2; k<c2+6; k++){ // defines the column
			for (int j=c2; j<c2+6; j++)// defines the row
				P.coeffRef(j-c2+6, k-c2+6) = s[k*ncols+j];
			for (int j=c1; j<c1+6; j++)
				P.coeffRef(j-c1, k-c2+6) = s[k*ncols+j];
		}; //std::cout << '\n' << P << std::endl;

		// measurement covariance
		for (int j=0; j<6; j++){
			for (int i=0; i<6; i++){
				R.coeffRef(i,j) = constraints[ii].R[j*6+i];
			}
		}; //std::cout << '\n' << R << std::endl;

		// jacobian
		constraint_jacobian_nodepair_model(Hs, xs, c1, c2);
		H = Eigen::MatrixXd(Hs);
		//std::cout << '\n' << H << std::endl;

		// predicted measurements
		std::vector<double> zs(6,0.0);
		constraint_model(zs, xs, c1, c2);
		//std::cout << '\n';
		//for (int i=0; i<6; i++)
		//  std::cout << zs[i] << std::endl;

		// innovations
		v = H*dx;
		for (int i=0; i<6; i++)
			v.coeffRef(i) = constraints[ii].z[i] - (zs[i] + v.coeffRef(i));
		pi_to_pi(v.coeffRef(4));
		//std::cout << '\n' << v << std::endl;

		// innovations covariance
		S = H*P*H.transpose() + R;
		//std::cout << '\n' << S << std::endl;

		// gate function
		gate[ii] = v.transpose()*(S.inverse()*v);
	}
}

/*=====================*/
void GraphOptimiser::constraint_model(std::vector<double>& zs, double *xs, int c1, int c2){
	std::vector<double> x1(6,0.0), x2(6,0.0);
	for (int i=0; i<6; i++){
		x1[i] = xs[c1+i];
		x2[i] = xs[c2+i];
	}
	std::vector<double> R1(9,0.0), R2(9,0.0), R1t(9,0.0), RR(9,0.0);
	compute_rotation_matrix(R1,x1);
	transpose(R1t,R1,3,3);
	for (int i=0; i<3; i++)
		zs[i] = R1t[i+0]*(x2[0]-x1[0])+R1t[i+3]*(x2[1]-x1[1])+R1t[i+6]*(x2[2]-x1[2]);
	compute_rotation_matrix(R2,x2);
	multiply(RR,R1t,R2,3,3,3);
	//display_matrix(R1t,3,3);
	//display_matrix(R2,3,3);
	//display_matrix(RR,3,3);
	compute_rotation_angles(RR,zs);
}

/*=====================*/
void GraphOptimiser::constraint_jacobian_nodepair_model(Eigen::SparseMatrix<double> &H, double *x, int c1, int c2){
	std::vector<double> x1(6,0.0), x2(6,0.0), hi(72);
	std::vector<T> tripletList;
	for (int i=0; i<6; i++){
		x1[i] = x[c1+i];
		x2[i] = x[c2+i];
	}
	compute_constraint_model_derivatives(hi,x1,x2);
	for (int j=0; j<12; j++){ // columns
		for (int i=0; i<6; i++){ // rows
			tripletList.push_back(T(i, j, hi[6*j+i]));
		}
	}
	H.setFromTriplets(tripletList.begin(), tripletList.end());
}

/*=====================*/
void GraphOptimiser::compute_constraint_model_derivatives(
		std::vector<double> &h,
		std::vector<double> &x1,
		std::vector<double> &x2){

	double dx = x2[0]-x1[0];
	double dy = x2[1]-x1[1];
	double dz = x2[2]-x1[2];

	// compute rotation matrix and its derivatives
	std::vector<double> dRa(9,0.0),dRb(9,0.0),dRc(9,0.0),R1(9,0.0),R2(9,0.0);
	compute_rotation_matrix_derivatives(R1,dRa,dRb,dRc,x1);
	//compute_rotation_matrix(R1,x1);
	compute_rotation_matrix(R2,x2);

	// xi
	double dxdxi=-R1[0], dxdyi=-R1[1], dxdzi=-R1[2];   // - transposed R
	double dydxi=-R1[3], dydyi=-R1[4], dydzi=-R1[5];
	double dzdxi=-R1[6], dzdyi=-R1[7], dzdzi=-R1[8];
	double dadxi=0, dadyi=0, dadzi=0;
	double dbdxi=0, dbdyi=0, dbdzi=0;
	double dcdxi=0, dcdyi=0, dcdzi=0;
	double dxdai=dRa[0]*dx+dRa[3]*dy+dRa[6]*dz;
	double dydai=dRa[1]*dx+dRa[4]*dy+dRa[7]*dz;
	double dzdai=dRa[2]*dx+dRa[5]*dy+dRa[8]*dz;
	double dxdbi=dRb[0]*dx+dRb[3]*dy+dRb[6]*dz;
	double dydbi=dRb[1]*dx+dRb[4]*dy+dRb[7]*dz;
	double dzdbi=dRb[2]*dx+dRb[5]*dy+dRb[8]*dz;
	double dxdci=dRc[0]*dx+dRc[3]*dy+dRc[6]*dz;
	double dydci=dRc[1]*dx+dRc[4]*dy+dRc[7]*dz;
	double dzdci=dRc[2]*dx+dRc[5]*dy+dRc[8]*dz;

	// xj
	double dxdxj=R1[0], dxdyj=R1[1], dxdzj=R1[2];  // transposed R
	double dydxj=R1[3], dydyj=R1[4], dydzj=R1[5];
	double dzdxj=R1[6], dzdyj=R1[7], dzdzj=R1[8];
	double dadxj=0, dadyj=0, dadzj=0;
	double dbdxj=0, dbdyj=0, dbdzj=0;
	double dcdxj=0, dcdyj=0, dcdzj=0;
	double dxdaj=0, dxdbj=0, dxdcj=0;
	double dydaj=0, dydbj=0, dydcj=0;
	double dzdaj=0, dzdbj=0, dzdcj=0;

	// Attitude terms
	double dadai=-1, dadbi=0, dadci=0;
	double dbdai=0, dbdbi=-1, dbdci=0;
	double dcdai=0, dcdbi=0, dcdci=-1;
	double dadaj=1, dadbj=0, dadcj=0;
	double dbdaj=0, dbdbj=1, dbdcj=0;
	double dcdaj=0, dcdbj=0, dcdcj=1;
	dx = x2[3]-x1[3];
	dy = x2[4]-x1[4];
	dz = x2[5]-x1[5];
	if (sqrt(dx*dx+dy*dy+dz*dz) > 0) { //norm(a)
		std::vector<double> dRRdaij(54),dRRdaijt(54),R1t(9),RR(9);
		std::vector<double> dadRR(27),dadRRt(27),dadaij(18);
		derivative_R1tR2(dRRdaij,x1,x2); //9x6
		//display_matrix(dRRdaij, 9, 6);
		transpose(R1t,R1,3,3);
		multiply(RR,R1t,R2,3,3,3);
		//display_matrix(RR,3,3);
		derivative_R2w(dadRR,RR); //9x3
		//display_matrix(dadRR,9,3);
		transpose(dadRRt,dadRR,9,3);
		multiply(dadaij,dadRRt,dRRdaij,3,9,6);
		//display_matrix(dadaij,3,6);
		///////////////////
		//I guess the symbolic derivation was with reversed inputs,
		//hence the extra minus sign
		////////////////////
		dadai=-dadaij[0];    dadbi=-dadaij[3];    dadci=-dadaij[6];
		dbdai=-dadaij[1];	dbdbi=-dadaij[4];    dbdci=-dadaij[7];
		dcdai=-dadaij[2];	dcdbi=-dadaij[5];    dcdci=-dadaij[8];
		dadaj=-dadaij[9];    dadbj=-dadaij[12];   dadcj=-dadaij[15];
		dbdaj=-dadaij[10];   dbdbj=-dadaij[13];   dbdcj=-dadaij[16];
		dcdaj=-dadaij[11];   dcdbj=-dadaij[14];   dcdcj=-dadaij[17];
	}
	h[0]=dxdxi;h[6]=dxdyi;h[12]=dxdzi;h[18]=dxdai;h[24]=dxdbi;h[30]=dxdci;h[36]=dxdxj;h[42]=dxdyj;h[48]=dxdzj;h[54]=dxdaj;h[60]=dxdbj;h[66]=dxdcj;
	h[1]=dydxi;h[7]=dydyi;h[13]=dydzi;h[19]=dydai;h[25]=dydbi;h[31]=dydci;h[37]=dydxj;h[43]=dydyj;h[49]=dydzj;h[55]=dydaj;h[61]=dydbj;h[67]=dydcj;
	h[2]=dzdxi;h[8]=dzdyi;h[14]=dzdzi;h[20]=dzdai;h[26]=dzdbi;h[32]=dzdci;h[38]=dzdxj;h[44]=dzdyj;h[50]=dzdzj;h[56]=dzdaj;h[62]=dzdbj;h[68]=dzdcj;
	h[3]=dadxi;h[9]=dadyi;h[15]=dadzi;h[21]=dadai;h[27]=dadbi;h[33]=dadci;h[39]=dadxj;h[45]=dadyj;h[51]=dadzj;h[57]=dadaj;h[63]=dadbj;h[69]=dadcj;
	h[4]=dbdxi;h[10]=dbdyi;h[16]=dbdzi;h[22]=dbdai;h[28]=dbdbi;h[34]=dbdci;h[40]=dbdxj;h[46]=dbdyj;h[52]=dbdzj;h[58]=dbdaj;h[64]=dbdbj;h[70]=dbdcj;
	h[5]=dcdxi;h[11]=dcdyi;h[17]=dcdzi;h[23]=dcdai;h[29]=dcdbi;h[35]=dcdci;h[41]=dcdxj;h[47]=dcdyj;h[53]=dcdzj;h[59]=dcdaj;h[65]=dcdbj;h[71]=dcdcj;
}

/*=====================*/
void GraphOptimiser::compute_rotation_angles(std::vector<double> &R, std::vector<double> &x){
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
void GraphOptimiser::compute_rotation_matrix(std::vector<double> &R, std::vector<double> &x){
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
void GraphOptimiser::compute_rotation_matrix_derivatives(
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
void GraphOptimiser::derivative_R1tR2(std::vector<double> &dRRdaij,
		std::vector<double> &x1,
		std::vector<double> &x2){

	double a1=x1[3],b1=x1[4],c1=x1[5];
	double a2=x2[3],b2=x2[4],c2=x2[5];
	double n1 = sqrt(a1*a1+b1*b1+c1*c1)+0.000001, n12 = n1*n1, n13 = n12*n1, n14 = n12*n12;
	double n2 = sqrt(a2*a2+b2*b2+c2*c2)+0.000001, n22 = n2*n2, n23 = n22*n2, n24 = n22*n22;
	double sn1 = sin(n1);
	double cn1 = cos(n1);
	double sn2 = sin(n2);
	double cn2 = cos(n2);
	// drotda1
	dRRdaij[0] = ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*a1*c1)/n13 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14) - (cn2 - (a2*a2*(cn2 - 1))/(n22))*((2*a1*(cn1 - 1))/n12 + (sn1*a1)/n1 - (2*a1*a1*a1*(cn1 - 1))/n14 - (sn1*a1*a1*a1)/n13);
	dRRdaij[1] = ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((2*a1*b1*b1*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(b1*b1))/n13) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(a1*a1))/n12 - (sn1*(a1*a1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) + (cn2 - (a2*a2*(cn2 - 1))/(n22))*((sn1*a1*c1)/n13 - (cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[2] = ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((sn1*(a1*a1))/n13 - (cn1*(a1*a1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((2*a1*(c1*c1)*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(c1*c1))/n13) + (cn2 - (a2*a2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[3] = ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((2*a1*(cn1 - 1))/n12 + (sn1*a1)/n1 - (2*a1*a1*a1*(cn1 - 1))/n14 - (sn1*a1*a1*a1)/n13) + ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14) + (cn2 - (b2*b2*(cn2 - 1))/(n22))*((cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*a1*c1)/n13 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[4] = (cn2 - (b2*b2*(cn2 - 1))/(n22))*((2*a1*(b1*b1)*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(b1*b1))/n13) + ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(a1*a1))/n12 - (sn1*(a1*a1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((sn1*a1*c1)/n13 - (cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[5] = ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((2*a1*(c1*c1)*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(c1*c1))/n13) + (cn2 - (b2*b2*(cn2 - 1))/(n22))*((sn1*(a1*a1))/n13 - (cn1*(a1*a1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[6] = (cn2 - (c2*c2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*a1*c1)/n13 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14) - ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((2*a1*(cn1 - 1))/n12 + (sn1*a1)/n1 - (2*a1*a1*a1*(cn1 - 1))/n14 - (sn1*a1*a1*a1)/n13);
	dRRdaij[7] = (cn2 - (c2*c2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(a1*a1))/n12 - (sn1*(a1*a1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((2*a1*(b1*b1)*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(b1*b1))/n13) + ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((sn1*a1*c1)/n13 - (cn1*a1*c1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*a1*a1*b1)/n13 + (2*a1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[8] = (cn2 - (c2*c2*(cn2 - 1))/(n22))*((2*a1*(c1*c1)*(cn1 - 1))/n14 - (sn1*a1)/n1 + (sn1*a1*(c1*c1))/n13) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((sn1*(a1*a1))/n13 - (cn1*(a1*a1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) + ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*a1*a1*c1)/n13 + (2*a1*a1*c1*(cn1 - 1))/n14);
	// drotdb1
	dRRdaij[9] = (cn2 - (a2*a2*(cn2 - 1))/(n22))*((2*b1*(a1*a1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(a1*a1))/n13) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((sn1*(b1*b1))/n13 - (cn1*(b1*b1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) + ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*b1*c1)/n13 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[10] = (cn2 - (a2*a2*(cn2 - 1))/(n22))*((sn1*b1*c1)/n13 - (cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14) - ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((2*b1*(cn1 - 1))/n12 + (sn1*b1)/n1 - (2*b1*b1*b1*(cn1 - 1))/n14 - (sn1*b1*b1*b1)/n13);
	dRRdaij[11] = (cn2 - (a2*a2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(b1*b1))/n12 - (sn1*(b1*b1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((2*b1*(c1*c1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(c1*c1))/n13) + ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[12] = ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((sn1*(b1*b1))/n13 - (cn1*(b1*b1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((2*b1*(a1*a1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(a1*a1))/n13) + (cn2 - (b2*b2*(cn2 - 1))/(n22))*((cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*b1*c1)/n13 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[13] = ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((sn1*b1*c1)/n13 - (cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14) - (cn2 - (b2*b2*(cn2 - 1))/(n22))*((2*b1*(cn1 - 1))/n12 + (sn1*b1)/n1 - (2*b1*b1*b1*(cn1 - 1))/n14 - (sn1*b1*b1*b1)/n13);
	dRRdaij[14] = ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((2*b1*(c1*c1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(c1*c1))/n13) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(b1*b1))/n12 - (sn1*(b1*b1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) + (cn2 - (b2*b2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[15] = ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((2*b1*(a1*a1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(a1*a1))/n13) + (cn2 - (c2*c2*(cn2 - 1))/(n22))*((sn1*(b1*b1))/n13 - (cn1*(b1*b1))/n12 - sn1/n1 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*b1*c1)/n13 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14);
	dRRdaij[16] = ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((2*b1*(cn1 - 1))/n12 + (sn1*b1)/n1 - (2*b1*b1*b1*(cn1 - 1))/n14 - (sn1*b1*b1*b1)/n13) + ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((sn1*b1*c1)/n13 - (cn1*b1*c1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*b1*a1*b1)/n13 + (2*b1*a1*b1*(cn1 - 1))/n14) + (cn2 - (c2*c2*(cn2 - 1))/(n22))*((cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 - (sn1*b1*a1)/n13 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[17] = (cn2 - (c2*c2*(cn2 - 1))/(n22))*((2*b1*(c1*c1)*(cn1 - 1))/n14 - (sn1*b1)/n1 + (sn1*b1*(c1*c1))/n13) + ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(b1*b1))/n12 - (sn1*(b1*b1))/n13 + (sn1*a1*b1*c1)/n13 + (2*a1*b1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((sn1*b1*a1)/n13 - (cn1*b1*a1)/n12 - (c1*(cn1 - 1))/n12 + (sn1*b1*b1*c1)/n13 + (2*b1*b1*c1*(cn1 - 1))/n14);
	// drotdc1
	dRRdaij[18] = (cn2 - (a2*a2*(cn2 - 1))/(n22))*((2*c1*(a1*a1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(a1*a1))/n13) + ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(c1*c1))/n12 - (sn1*(c1*c1))/n13 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((sn1*c1*b1)/n13 - (cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[19] = ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((2*c1*(b1*b1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(b1*b1))/n13) + (cn2 - (a2*a2*(cn2 - 1))/(n22))*((sn1*(c1*c1))/n13 - (cn1*(c1*c1))/n12 - sn1/n1 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) - ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*c1*a1)/n13 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[20] = ((b2*sn2)/n2 + (a2*c2*(cn2 - 1))/(n22))*((2*c1*(cn1 - 1))/n12 + (sn1*c1)/n1 - (2*c1*c1*c1*(cn1 - 1))/n14 - (sn1*c1*c1*c1)/n13) + ((c2*sn2)/n2 - (a2*b2*(cn2 - 1))/(n22))*((sn1*c1*a1)/n13 - (cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14) + (cn2 - (a2*a2*(cn2 - 1))/(n22))*((cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*c1*b1)/n13 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[21] = (cn2 - (b2*b2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(c1*c1))/n12 - (sn1*(c1*c1))/n13 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((2*c1*(a1*a1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(a1*a1))/n13) + ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((sn1*c1*b1)/n13 - (cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[22] = (cn2 - (b2*b2*(cn2 - 1))/(n22))*((2*c1*(b1*b1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(b1*b1))/n13) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((sn1*(c1*c1))/n13 - (cn1*(c1*c1))/n12 - sn1/n1 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) + ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*c1*a1)/n13 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[23] = (cn2 - (b2*b2*(cn2 - 1))/(n22))*((sn1*c1*a1)/n13 - (cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 - (b2*c2*(cn2 - 1))/(n22))*((2*c1*(cn1 - 1))/n12 + (sn1*c1)/n1 - (2*c1*c1*c1*(cn1 - 1))/n14 - (sn1*c1*c1*c1)/n13) - ((c2*sn2)/n2 + (a2*b2*(cn2 - 1))/(n22))*((cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*c1*b1)/n13 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[24] = ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((2*c1*(a1*a1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(a1*a1))/n13) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*(sn1/n1 + (cn1*(c1*c1))/n12 - (sn1*(c1*c1))/n13 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) + (cn2 - (c2*c2*(cn2 - 1))/(n22))*((sn1*c1*b1)/n13 - (cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14);
	dRRdaij[25] = ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((sn1*(c1*c1))/n13 - (cn1*(c1*c1))/n12 - sn1/n1 + (sn1*c1*a1*b1)/n13 + (2*c1*a1*b1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((2*c1*(b1*b1)*(cn1 - 1))/n14 - (sn1*c1)/n1 + (sn1*c1*(b1*b1))/n13) + (cn2 - (c2*c2*(cn2 - 1))/(n22))*((cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 - (sn1*c1*a1)/n13 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14);
	dRRdaij[26] = ((b2*sn2)/n2 - (a2*c2*(cn2 - 1))/(n22))*((cn1*c1*b1)/n12 - (a1*(cn1 - 1))/n12 - (sn1*c1*b1)/n13 + (sn1*c1*a1*c1)/n13 + (2*c1*a1*c1*(cn1 - 1))/n14) - ((a2*sn2)/n2 + (b2*c2*(cn2 - 1))/(n22))*((sn1*c1*a1)/n13 - (cn1*c1*a1)/n12 - (b1*(cn1 - 1))/n12 + (sn1*c1*b1*c1)/n13 + (2*c1*b1*c1*(cn1 - 1))/n14) - (cn2 - (c2*c2*(cn2 - 1))/(n22))*((2*c1*(cn1 - 1))/n12 + (sn1*c1)/n1 - (2*c1*c1*c1*(cn1 - 1))/n14 - (sn1*c1*c1*c1)/n13);
	// drotda2
	dRRdaij[27] = ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) - (c2*sn2*a2)/n23 + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) - (cn1 - (a1*a1*(cn1 - 1))/n12)*((2*a2*(cn2 - 1))/(n22) + (sn2*a2)/n2 - (a2*a2*sn2*a2)/n23 - (2*a2*a2*a2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((b2*sn2*a2)/n23 - (b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24);
	dRRdaij[28] = ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((b2*sn2*a2)/n23 - (b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24) + (cn1 - (b1*b1*(cn1 - 1))/n12)*((c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) - (c2*sn2*a2)/n23 + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) + ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((2*a2*(cn2 - 1))/(n22) + (sn2*a2)/n2 - (a2*a2*sn2*a2)/n23 - (2*a2*a2*a2*(cn2 - 1))/n24);
	dRRdaij[29] = (cn1 - (c1*c1*(cn1 - 1))/n12)*((b2*sn2*a2)/n23 - (b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) - (c2*sn2*a2)/n23 + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) - ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((2*a2*(cn2 - 1))/(n22) + (sn2*a2)/n2 - (a2*a2*sn2*a2)/n23 - (2*a2*a2*a2*(cn2 - 1))/n24);
	dRRdaij[30] = (cn1 - (a1*a1*(cn1 - 1))/n12)*((c2*sn2*a2)/n23 - (c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*(sn2/n2 + (a2*cn2*a2)/(n22) - (a2*sn2*a2)/n23 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24) + ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((b2*b2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*b2*b2*a2*(cn2 - 1))/n24);
	dRRdaij[31] = (cn1 - (b1*b1*(cn1 - 1))/n12)*((b2*b2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*b2*b2*a2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((c2*sn2*a2)/n23 - (c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) + ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*(sn2/n2 + (a2*cn2*a2)/(n22) - (a2*sn2*a2)/n23 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24);
	dRRdaij[32] = ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((c2*sn2*a2)/n23 - (c2*cn2*a2)/(n22) - (b2*(cn2 - 1))/(n22) + (a2*b2*sn2*a2)/n23 + (2*a2*b2*a2*(cn2 - 1))/n24) + (cn1 - (c1*c1*(cn1 - 1))/n12)*(sn2/n2 + (a2*cn2*a2)/(n22) - (a2*sn2*a2)/n23 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((b2*b2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*b2*b2*a2*(cn2 - 1))/n24);
	dRRdaij[33] = (cn1 - (a1*a1*(cn1 - 1))/n12)*((b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) - (b2*sn2*a2)/n23 + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24) + ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((a2*sn2*a2)/n23 - (a2*cn2*a2)/(n22) - sn2/n2 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*c2*c2*a2*(cn2 - 1))/n24);
	dRRdaij[34] = (cn1 - (b1*b1*(cn1 - 1))/n12)*((a2*sn2*a2)/n23 - (a2*cn2*a2)/(n22) - sn2/n2 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) - (b2*sn2*a2)/n23 + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24) + ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*c2*c2*a2*(cn2 - 1))/n24);
	dRRdaij[35] = (cn1 - (c1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*a2)/n23 - (sn2*a2)/n2 + (2*c2*c2*a2*(cn2 - 1))/n24) + ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((b2*cn2*a2)/(n22) - (c2*(cn2 - 1))/(n22) - (b2*sn2*a2)/n23 + (a2*c2*sn2*a2)/n23 + (2*a2*c2*a2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((a2*sn2*a2)/n23 - (a2*cn2*a2)/(n22) - sn2/n2 + (b2*c2*sn2*a2)/n23 + (2*b2*c2*a2*(cn2 - 1))/n24);
	// drotdb2
	dRRdaij[36] = (cn1 - (a1*a1*(cn1 - 1))/n12)*((a2*a2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*a2*a2*b2*(cn2 - 1))/n24) + ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) - (c2*sn2*b2)/n23 + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((b2*sn2*b2)/n23 - (b2*cn2*b2)/(n22) - sn2/n2 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	dRRdaij[37] = (cn1 - (b1*b1*(cn1 - 1))/n12)*((c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) - (c2*sn2*b2)/n23 + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((a2*a2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*a2*a2*b2*(cn2 - 1))/n24) + ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((b2*sn2*b2)/n23 - (b2*cn2*b2)/(n22) - sn2/n2 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	dRRdaij[38] = ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((a2*a2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*a2*a2*b2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) - (c2*sn2*b2)/n23 + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) + (cn1 - (c1*c1*(cn1 - 1))/n12)*((b2*sn2*b2)/n23 - (b2*cn2*b2)/(n22) - sn2/n2 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	dRRdaij[39] = (cn1 - (a1*a1*(cn1 - 1))/n12)*((c2*sn2*b2)/n23 - (c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) - (a2*sn2*b2)/n23 + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) - ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((2*b2*(cn2 - 1))/(n22) + (sn2*b2)/n2 - (b2*b2*sn2*b2)/n23 - (2*b2*b2*b2*(cn2 - 1))/n24);
	dRRdaij[40] = ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) - (a2*sn2*b2)/n23 + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((c2*sn2*b2)/n23 - (c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) - (cn1 - (b1*b1*(cn1 - 1))/n12)*((2*b2*(cn2 - 1))/(n22) + (sn2*b2)/n2 - (b2*b2*sn2*b2)/n23 - (2*b2*b2*b2*(cn2 - 1))/n24);
	dRRdaij[41] = ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((c2*sn2*b2)/n23 - (c2*cn2*b2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*b2*sn2*b2)/n23 + (2*a2*b2*b2*(cn2 - 1))/n24) + (cn1 - (c1*c1*(cn1 - 1))/n12)*((a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) - (a2*sn2*b2)/n23 + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) + ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((2*b2*(cn2 - 1))/(n22) + (sn2*b2)/n2 - (b2*b2*sn2*b2)/n23 - (2*b2*b2*b2*(cn2 - 1))/n24);
	dRRdaij[42] = ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((a2*sn2*b2)/n23 - (a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*c2*c2*b2*(cn2 - 1))/n24) + (cn1 - (a1*a1*(cn1 - 1))/n12)*(sn2/n2 + (b2*cn2*b2)/(n22) - (b2*sn2*b2)/n23 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	dRRdaij[43] = ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*c2*c2*b2*(cn2 - 1))/n24) + (cn1 - (b1*b1*(cn1 - 1))/n12)*((a2*sn2*b2)/n23 - (a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*(sn2/n2 + (b2*cn2*b2)/(n22) - (b2*sn2*b2)/n23 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	dRRdaij[44] = (cn1 - (c1*c1*(cn1 - 1))/n12)*((c2*c2*sn2*b2)/n23 - (sn2*b2)/n2 + (2*c2*c2*b2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((a2*sn2*b2)/n23 - (a2*cn2*b2)/(n22) - (c2*(cn2 - 1))/(n22) + (b2*c2*sn2*b2)/n23 + (2*b2*c2*b2*(cn2 - 1))/n24) + ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*(sn2/n2 + (b2*cn2*b2)/(n22) - (b2*sn2*b2)/n23 + (a2*c2*sn2*b2)/n23 + (2*a2*c2*b2*(cn2 - 1))/n24);
	// drotdc2
	dRRdaij[45] = (cn1 - (a1*a1*(cn1 - 1))/n12)*((a2*a2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*a2*a2*c2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((b2*sn2*c2)/n23 - (b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) + ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*(sn2/n2 + (c2*cn2*c2)/(n22) - (c2*sn2*c2)/n23 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[46] = ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((b2*sn2*c2)/n23 - (b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((a2*a2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*a2*a2*c2*(cn2 - 1))/n24) + (cn1 - (b1*b1*(cn1 - 1))/n12)*(sn2/n2 + (c2*cn2*c2)/(n22) - (c2*sn2*c2)/n23 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[47] = ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((a2*a2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*a2*a2*c2*(cn2 - 1))/n24) + (cn1 - (c1*c1*(cn1 - 1))/n12)*((b2*sn2*c2)/n23 - (b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*(sn2/n2 + (c2*cn2*c2)/(n22) - (c2*sn2*c2)/n23 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[48] = ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((b2*b2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*b2*b2*c2*(cn2 - 1))/n24) - ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) - (a2*sn2*c2)/n23 + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24) + (cn1 - (a1*a1*(cn1 - 1))/n12)*((c2*sn2*c2)/n23 - (c2*cn2*c2)/(n22) - sn2/n2 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[49] = (cn1 - (b1*b1*(cn1 - 1))/n12)*((b2*b2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*b2*b2*c2*(cn2 - 1))/n24) + ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) - (a2*sn2*c2)/n23 + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((c2*sn2*c2)/n23 - (c2*cn2*c2)/(n22) - sn2/n2 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[50] = (cn1 - (c1*c1*(cn1 - 1))/n12)*((a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) - (a2*sn2*c2)/n23 + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((b2*b2*sn2*c2)/n23 - (sn2*c2)/n2 + (2*b2*b2*c2*(cn2 - 1))/n24) + ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((c2*sn2*c2)/n23 - (c2*cn2*c2)/(n22) - sn2/n2 + (a2*b2*sn2*c2)/n23 + (2*a2*b2*c2*(cn2 - 1))/n24);
	dRRdaij[51] = ((sn1*c1)/n1 - (a1*b1*(cn1 - 1))/n12)*((a2*sn2*c2)/n23 - (a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24) + (cn1 - (a1*a1*(cn1 - 1))/n12)*((b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) - (b2*sn2*c2)/n23 + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) + ((sn1*b1)/n1 + (a1*c1*(cn1 - 1))/n12)*((2*c2*(cn2 - 1))/(n22) + (sn2*c2)/n2 - (c2*c2*sn2*c2)/n23 - (2*c2*c2*c2*(cn2 - 1))/n24);
	dRRdaij[52] = (cn1 - (b1*b1*(cn1 - 1))/n12)*((a2*sn2*c2)/n23 - (a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24) - ((sn1*c1)/n1 + (a1*b1*(cn1 - 1))/n12)*((b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) - (b2*sn2*c2)/n23 + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) - ((sn1*a1)/n1 - (b1*c1*(cn1 - 1))/n12)*((2*c2*(cn2 - 1))/(n22) + (sn2*c2)/n2 - (c2*c2*sn2*c2)/n23 - (2*c2*c2*c2*(cn2 - 1))/n24);
	dRRdaij[53] = ((sn1*b1)/n1 - (a1*c1*(cn1 - 1))/n12)*((b2*cn2*c2)/(n22) - (a2*(cn2 - 1))/(n22) - (b2*sn2*c2)/n23 + (a2*c2*sn2*c2)/n23 + (2*a2*c2*c2*(cn2 - 1))/n24) - (cn1 - (c1*c1*(cn1 - 1))/n12)*((2*c2*(cn2 - 1))/(n22) + (sn2*c2)/n2 - (c2*c2*sn2*c2)/n23 - (2*c2*c2*c2*(cn2 - 1))/n24) - ((sn1*a1)/n1 + (b1*c1*(cn1 - 1))/n12)*((a2*sn2*c2)/n23 - (a2*cn2*c2)/(n22) - (b2*(cn2 - 1))/(n22) + (b2*c2*sn2*c2)/n23 + (2*b2*c2*c2*(cn2 - 1))/n24);
}


/*=====================*/
void GraphOptimiser::derivative_R2w(std::vector<double> &dadRR,
		std::vector<double> &R){
	double r1=R[0],r2=R[1],r3=R[2];
	double r4=R[3],r5=R[4],r6=R[5];
	double r7=R[6],r8=R[7],r9=R[8];
	double a = ((r2/2-r4/2)*(r2/2-r4/2)+(r3/2-r7/2)*(r3/2-r7/2)+(r6/2-r8/2)*(r6/2-r8/2));
	double a2= sqrt(a*a*a);
	double b = 0.5*(r1+r5+r9-1);
	double c = (r6/2-r8/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
	double d = (r3/2-r7/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
	double e = (r2/2-r4/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
	double f = (r2/2-r4/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
	double g = (r2/2-r4/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
	double h = (r2/2-r4/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
	double i = (r2/2-r4/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
	double j = atan2(sqrt(a),b) / (2*sqrt(a));
	double k = (r2/2-r4/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
	double l = (r2/2-r4/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
	double m = (r3/2-r7/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
	double n = (r3/2-r7/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
	double o = (r3/2-r7/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
	double p = (r3/2-r7/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
	double q = (r3/2-r7/2)*atan2(sqrt(a),b)*(r2/2-r4/2)/(2*a2);
	double r = (r3/2-r7/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
	double s = (r6/2-r8/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
	double t = (r6/2-r8/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
	double u = (r6/2-r8/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
	double v = (r6/2-r8/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
	double w = (r6/2-r8/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
	double x = (r6/2-r8/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
	double y = (r3/2-r7/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
	double z =(r6/2 -r8/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
	dadRR[0] = c;       dadRR[9] = -d;      dadRR[18] = e;
	dadRR[1] = f-g;     dadRR[10] = h-i;	dadRR[19] = -j+k-l;
	dadRR[2] = m-n;     dadRR[11] = j-o+p;  dadRR[20] = q-r;
	dadRR[3] = g-f;     dadRR[12] = i- h;	dadRR[21] = j-k+l;
	dadRR[4] = c;       dadRR[13] = -d;     dadRR[22] = e;
	dadRR[5] = -j+s-t;	dadRR[14] = u-v;	dadRR[23] = w-x;
	dadRR[6] = n-m;     dadRR[15] = -j+o-p; dadRR[24] = r-y;
	dadRR[7] = j-s+t;   dadRR[16] = v-z;	dadRR[25] = x-w;
	dadRR[8] = c;       dadRR[17] = -d;     dadRR[26] = e;
}

/*=====================*/
void GraphOptimiser::multiply(std::vector<double> &C, std::vector<double> &A, std::vector<double> &B, int nrows, int ncomm, int ncols){
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
void GraphOptimiser::transpose(std::vector<double> &At, std::vector<double> &A,
		int nrows, int ncols){
	for (int j=0; j<ncols; j++){
		for (int i=0; i<nrows; i++){
			At[i*ncols+j] = A[ j*nrows+i];
		}
	}
}

/*=====================*/
void GraphOptimiser::display_matrix(std::vector<double> &A, int nrows, int ncols){
	std::cout << '\n';
	for (int j=0; j<nrows; j++){
		for (int i=0; i<ncols; i++){
			std::cout << A[i*nrows+j] << ' ';
		}
		std::cout << '\n';
	}
}

/*=====================*/
void GraphOptimiser::pi_to_pi(double &angle){
	double pi = 3.1415926534;
	double twopi = 2*pi;
	//angle = angle - twopi*fix(angle/twopi); // this is a stripped-down version of rem(angle, 2*pi)
	angle = angle - twopi*(angle - remainder(angle,twopi))/twopi;
	if (angle > pi){
		angle = angle - twopi;};
	if (angle < -pi){
		angle = angle + twopi;};
}

/*
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */
