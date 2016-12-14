function H=relpose_jacobian_6DoF(x1,x2)
%function H=relpose_jacobian_6DoF(x)
%
% INPUT:  x - state vector of two poses
% OUTPUT: H - relative pose Jacobian matrix (about point x)
%
% Given two poses x(1:6) and x(7:12), the relative pose is given by
% xr = transform_to_relative(x(4:6), x(1:3)).
% Here we compute the Jacobian of that transformation.
%
% test data 
% note, code is slow at first run, due to compile after any code changes
% x=rand(12,1);
% tic; analytical=relpose_jacobian_6DoF(x(1:6), x(7:12)); toc
% tic; numerical=numerical_jacobian_i(@constraint_model, [], 1, [], x); toc
% error = analytical - numerical
%
% By Tariq Abuhashim 22/03/2016
%

dx=x2(1)-x1(1);
dy=x2(2)-x1(2);
dz=x2(3)-x1(3);
[R,dRa,dRb,dRc]=get_rotation_matrix(x1(4:6),'rodrigues');
R=R';
dRa=dRa';
dRb=dRb';
dRc=dRc';

% xi
dxdxi=-R(1,1);  dxdyi=-R(1,2);  dxdzi=-R(1,3);
dydxi=-R(2,1);  dydyi=-R(2,2);  dydzi=-R(2,3);
dzdxi=-R(3,1);  dzdyi=-R(3,2);  dzdzi=-R(3,3);
dadxi=0;	dadyi=0;	dadzi=0;
dbdxi=0;	dbdyi=0;	dbdzi=0;
dcdxi=0;	dcdyi=0;	dcdzi=0;
dxdai=dRa(1,:)*[dx;dy;dz];  dxdbi=dRb(1,:)*[dx;dy;dz];  dxdci=dRc(1,:)*[dx;dy;dz];
dydai=dRa(2,:)*[dx;dy;dz];  dydbi=dRb(2,:)*[dx;dy;dz];  dydci=dRc(2,:)*[dx;dy;dz];
dzdai=dRa(3,:)*[dx;dy;dz];  dzdbi=dRb(3,:)*[dx;dy;dz];  dzdci=dRc(3,:)*[dx;dy;dz];

% xj
dxdxj=R(1,1);   dxdyj=R(1,2);   dxdzj=R(1,3);
dydxj=R(2,1);   dydyj=R(2,2);   dydzj=R(2,3);
dzdxj=R(3,1);   dzdyj=R(3,2);   dzdzj=R(3,3);
dadxj=0;	dadyj=0;	dadzj=0;
dbdxj=0;	dbdyj=0;	dbdzj=0;
dcdxj=0;	dcdyj=0;	dcdzj=0;
dxdaj=0;    dxdbj=0;    dxdcj=0;
dydaj=0;    dydbj=0;    dydcj=0;
dzdaj=0;    dzdbj=0;    dzdcj=0;

%%%%%%%% Attitude terms
% check the file analytical_derivation_of_rotation_jacobians.m
% for the use of MATLAB Symbolic toolbox for derivation
R1=w2R(x1(4:6));
R2=w2R(x2(4:6));
RR=R1'*R2;
a=R2w(RR);
a(2)=pi_to_pi(a(2));
if norm(a)
    %dRRdaij = derivative_R1tR2_mex(x1(4:6),x2(4:6)); %9x6
    dRRdaij = derivative_R1tR2(x1(4:6),x2(4:6)); %9x6
    dadRR = derivative_R2w(RR); %9x3
    dadaij = -dadRR'*dRRdaij; %3x6
    %I guess the symbolic derivation was with reversed inputs,
    %hence the extra minus sign
    dadai=dadaij(1);	dadbi=dadaij(4);	dadci=dadaij(7);
    dbdai=dadaij(2);	dbdbi=dadaij(5);	dbdci=dadaij(8);
    dcdai=dadaij(3);	dcdbi=dadaij(6);	dcdci=dadaij(9);
    dadaj=dadaij(10);    dadbj=dadaij(13);    dadcj=dadaij(16);
    dbdaj=dadaij(11);    dbdbj=dadaij(14);    dbdcj=dadaij(17);
    dcdaj=dadaij(12);    dcdbj=dadaij(15);    dcdcj=dadaij(18);
else
    dadai=-1;	dadbi=0;	dadci=0;
    dbdai=0;	dbdbi=-1;	dbdci=0;
    dcdai=0;	dcdbi=0;	dcdci=-1;
    dadaj=1;    dadbj=0;    dadcj=0;
    dbdaj=0;    dbdbj=1;    dbdcj=0;
    dcdaj=0;    dcdbj=0;    dcdcj=1;
end
%%%%%%%%

H=[
    dxdxi dxdyi dxdzi,  dxdai dxdbi dxdci,  dxdxj dxdyj dxdzj,  dxdaj dxdbj dxdcj;
    dydxi dydyi dydzi,  dydai dydbi dydci,  dydxj dydyj dydzj,  dydaj dydbj dydcj;
    dzdxi dzdyi dzdzi,  dzdai dzdbi dzdci,  dzdxj dzdyj dzdzj,  dzdaj dzdbj dzdcj;
    
    dadxi dadyi dadzi,  dadai dadbi dadci,  dadxj dadyj dadzj,  dadaj dadbj dadcj;
    dbdxi dbdyi dbdzi,  dbdai dbdbi dbdci,  dbdxj dbdyj dbdzj,  dbdaj dbdbj dbdcj;
    dcdxi dcdyi dcdzi,  dcdai dcdbi dcdci,  dcdxj dcdyj dcdzj,  dcdaj dcdbj dcdcj
    ];

% % by using numerical jacobian?
% % transform_to_relative_w takes x2 as input 1 and x1 as input 2;
% J1 = numerical_jacobian_i(@transform_to_relative_w,[],2,[],x2,x1);
% J2 = numerical_jacobian_i(@transform_to_relative_w,[],1,[],x2,x1);
% %H = [J1 J2];
% sss = abs(H - [J1 J2]);
% if max(sss(:)) > 1e-10
%     [J1 J2]
%     H
%     pause
% end

function df=derivative_R2w(R)
%
% derivative were produced using:
%
% syms n(r1,r2,r3,r4,r5,r6,r7,r8,r9)
% n = R2w([r1,r2,r3;r4,r5,r6;r7,r8,r9]);
% diff(n,r1)
% diff(n,r2)
% diff(n,r3)
% diff(n,r4)
% diff(n,r5)
% diff(n,r6)
% diff(n,r7)
% diff(n,r8)
% diff(n,r9)
%
% By Tariq Abuhashim 22/03/2016
%

r1=R(1);r2=R(2);r3=R(3);
r4=R(4);r5=R(5);r6=R(6);
r7=R(7);r8=R(8);r9=R(9);
a = ((r2/2-r4/2)*(r2/2-r4/2)+(r3/2-r7/2)*(r3/2-r7/2)+(r6/2-r8/2)*(r6/2-r8/2));
a2= a^(3/2);
b = r1/2+r5/2+r9/2-1/2;
c = (r6/2-r8/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
d = (r3/2-r7/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
e = (r2/2-r4/2)*(sqrt(a)/2) / ((a+b*b)*sqrt(a));
f = (r2/2-r4/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
g = (r2/2-r4/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
h = (r2/2-r4/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
i = (r2/2-r4/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
j = atan2(sqrt(a),b) / (2*sqrt(a));
k = (r2/2-r4/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
l = (r2/2-r4/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
m = (r3/2-r7/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
n = (r3/2-r7/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
o = (r3/2-r7/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
p = (r3/2-r7/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
q = (r3/2-r7/2)*atan2(sqrt(a),b)*(r2/2-r4/2)/(2*a2);
r = (r3/2-r7/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
s = (r6/2-r8/2)*atan2(sqrt(a),b)*(r6/2-r8/2) / (2*a2);
t = (r6/2-r8/2)*(r6/2-r8/2)*b / (2*(a+b*b)*a);
u = (r6/2-r8/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);
v = (r6/2-r8/2)*atan2(sqrt(a),b)*(r3/2-r7/2) / (2*a2);
w = (r6/2-r8/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
x = (r6/2-r8/2)*(r2/2-r4/2)*b / (2*(a+b*b)*a);
y = (r3/2-r7/2)*atan2(sqrt(a),b)*(r2/2-r4/2) / (2*a2);
z =(r6/2 -r8/2)*(r3/2-r7/2)*b / (2*(a+b*b)*a);

df = ...
[c,...
-d,...
 e;
 f - g, ...
 h - i, ...
-j + k - l;
 m - n,...
 j - o + p, ...
 q - r;
 g - f, ...
 i - h, ...
 j - k + l;
 c, ...
-d, ...
 e;
-j + s - t, ...
 u - v, ...
 w - x;
 n - m, ...
-j + o - p, ...
 r - y;
 j - s + t, ...
 v - z, ...
 x - w;
 c, ...
-d, ...
 e];