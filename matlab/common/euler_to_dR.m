% DCM derivatives

function [dRr, dRp, dRq] = euler_to_dR(a)

% calculates from euler angles (a) the derivatives of the direction cosine 
% matrix (dRr, dRp, dRq).
% a: should be in randians, and in an order of [roll pitch yaw], 'XYZ'
% R = [a1 a2 a3]
%     [b1 b2 b3]
%     [c1 c2 c3]
%
% a = [roll pitch yaw] = [r p q]

c=cos(a);
s=sin(a);

% wrt roll angle
da1r = 0; 
da2r = c(3)*s(2)*c(1) + s(3)*s(1); 
da3r =-c(3)*s(2)*s(1) + s(3)*c(1);
db1r = 0; 
db2r = s(3)*s(2)*c(1) - c(3)*s(1); 
db3r =-s(3)*s(2)*s(1) - c(3)*c(1);
dc1r = 0; 
dc2r = c(2)*c(1);                
dc3r =-c(2)*s(1);

% wrt pitch angle
da1p =-c(3)*s(2); 
da2p = c(3)*c(2)*s(1); 
da3p = c(3)*c(2)*c(1);
db1p =-s(3)*s(2); 
db2p = s(3)*c(2)*s(1); 
db3p = s(3)*c(2)*c(1);
dc1p =-c(2);      
dc2p =-s(2)*s(1);      
dc3p =-s(2)*c(1);

% wrt roll angle
da1q =-s(3)*c(2);
da2q =-s(3)*s(2)*s(1) - c(3)*c(1);
da3q =-s(3)*s(2)*c(1) + c(3)*s(1);
db1q = c(3)*c(2);
db2q = c(3)*s(2)*s(1) - s(3)*c(1);
db3q = c(3)*s(2)*c(1) + s(3)*s(1);
dc1q = 0;
dc2q = 0;
dc3q = 0;

% all together
dRr = [da1r da2r da3r; db1r db2r db3r; dc1r dc2r dc3r]; % roll derivatives
dRp = [da1p da2p da3p; db1p db2p db3p; dc1p dc2p dc3p]; % pitch derivatives
dRq = [da1q da2q da3q; db1q db2q db3q; dc1q dc2q dc3q]; % yaw derivatives