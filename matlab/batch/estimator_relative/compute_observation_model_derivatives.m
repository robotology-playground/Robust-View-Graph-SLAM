function hi = compute_observation_model_derivatives(x, p, varargin)
if size(p,1) < 3; p = pextend(p); end;
% get rotation and its derivatives ('rodrigues', 'exp_skew', 'dcm')
[R, dRa, dRb, dRc] = get_rotation_matrix(x(4:6), 'rodrigues');
R = R';
dRa = dRa';
dRb = dRb';
dRc = dRc';

% range
r = 1./x(7); % range from inverse range

% 3d points in camera 1 frame
%[p, dp] = get_scan_from_range(p, r);
if size(p,1) < 3; p = pextend(p); end;
rim = sqrt(p(1)*p(1) + p(2)*p(2) + 1);
d = r./rim;
dd = -r.*d;
dp = dd.*p;
p = d*p;
% Derivatives of camera 2 measurements
%xc = x(end-5:end); % relative pose
%xr = transform_to_relative_w(xf,xc);
p = p - x(1:3);
num1 = R(1,:)*p;
num2 = R(2,:)*p;
denum1 = R(3,:)*p;
denum2 = denum1.^2;
dnum1 = R(1,:)*dp;
dnum2 = R(2,:)*dp;
ddenum1 = R(3,:)*dp;
%hi = zeros(14, 1);
% derivatives of pix with respect to inverse depth rho = x(1)
hi(1) = (dnum1.*denum1 - num1.*ddenum1)./denum2;
hi(2) = (dnum2.*denum1 - num2.*ddenum1)./denum2;
%hi(1) = R(1,:)*dxf./denum1 - num1*R(3,:)*dxf./denum2;
%hi(2) = R(2,:)*dxf./denum1 - num2*R(3,:)*dxf./denum2;
% derivatives of pix with respect to translation t = [t1 t2 t3]'.
temp = R(3,1)./denum2;
hi(3) = - R(1,1)./denum1 + num1.*temp; % dudt1
hi(4) = - R(2,1)./denum1 + num2.*temp; % dvdt1
temp = R(3,2)./denum2;
hi(5) = - R(1,2)./denum1 + num1.*temp; % dudt2
hi(6) = - R(2,2)./denum1 + num2.*temp; % dvdt2
temp = R(3,3)./denum2;
hi(7) = - R(1,3)./denum1 + num1.*temp; % dudt3
hi(8) = - R(2,3)./denum1 + num2.*temp; % dvdt3
% derivatives of pix with respect to rotation a = [a1 a2 a3]'.
xfa = dRa*p;
temp = xfa(3,:)./denum2;
hi(9)  = xfa(1,:)./denum1 - num1.*temp; % duda
hi(10) = xfa(2,:)./denum1 - num2.*temp; % dvda
xfb = dRb*p;
temp = xfb(3,:)./denum2;
hi(11) = xfb(1,:)./denum1 - num1.*temp; % dudb
hi(12) = xfb(2,:)./denum1 - num2.*temp; % dvdb
xfc = dRc*p;
temp = xfc(3,:)./denum2;
hi(13) = xfc(1,:)./denum1 - num1.*temp; % dudc
hi(14) = xfc(2,:)./denum1 - num2.*temp; % dvdc