#include "graph.h"
#include "averaging.h"

vector<double> averaging::RtoQuaternion(vector<double> R){
    Q(1)=(R.coeffRef(1,1)+R.coeffRef(2,2)+R.coeffRef(3,3)-1)/2;
    Q(2)=(R.coeffRef(3,2)-R.coeffRef(2,3))/2;
    Q(3)=(R.coeffRef(1,3)-R.coeffRef(3,1))/2;
    Q(4)=(R.coeffRef(2,1)-R.coeffRef(1,2))/2;
    Q(1)=sqrt((Q(1)+1)/2);
    Q(2)=(Q(2)/Q(1))/2;
    Q(3)=(Q(3)/Q(1))/2;
    Q(4)=(Q(4)/Q(1))/2;}

void averaging::QuaternionAvg(vector<Graph::Constraints> C){
    std::cout << C.size() << endl;
}

vector<double> averaging::QuaternionFwd(vector<double> Qi, vector<double> Qij){
    // Qj=Qij*Qi
    vector<double> Qj[4];
    Qj[0]=Qij[0]*Qi[0]-Qij[1]*Qi[1]-Qij[2]*Qi[2]-Qij[3]*Qi[3];
    Qj[1]=Qij[1]*Qi[0]+Qij[0]*Qi[1]-Qij[3]*Qi[2]+Qij[2]*Qi[3];
    Qj[2]=Qij[2]*Qi[0]+Qij[3]*Qi[1]+Qij[0]*Qi[2]-Qij[1]*Qi[3];
    Qj[3]=Qij[3]*Qi[0]-Qij[2]*Qi[1]+Qij[1]*Qi[2]+Qij[0]*Qi[3];}

vector<double> averaging::QuaternionRvs(vector<double> Qj, vector<double> Qij){
    // Qi=Qij*Qj
    vector<double> Qi[4];
    Qi[0]=-Qij[0]*Qj[0]-Qij[1]*Qj[1]-Qij[2]*Qj[2]-Qij[3]*Qj[3];
    Qi[1]= Qij[1]*Qj[0]-Qij[0]*Qj[1]-Qij[3]*Qj[2]+Qij[2]*Qj[3];
    Qi[2]= Qij[2]*Qj[0]+Qij[3]*Qj[1]-Qij[0]*Qj[2]-Qij[1]*Qj[3];
    Qi[3]= Qij[3]*Qj[0]-Qij[2]*Qj[1]+Qij[1]*Qj[2]-Qij[0]*Qj[3];}