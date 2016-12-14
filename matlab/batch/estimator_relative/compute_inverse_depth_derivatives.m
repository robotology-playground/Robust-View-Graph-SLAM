function hi = compute_inverse_depth_derivatives(x1, x2, xs, varargin)
if size(x1,1) < 3; x1 = pextend(x1); end
if size(x2,1) < 3; x2 = pextend(x2); end
n = size(x1,2);
% get rotation and its derivatives ('rodrigues', 'exp_skew', 'dcm')
if isempty(varargin)
    [R, dR1, dR2, dR3] = get_rotation_matrix(xs(4:6), 'rodrigues');
else
	R = varargin{1};
	dR1 = varargin{2};
	dR2 = varargin{3};
	dR3 = varargin{4};
end
% calculate the terms
v1 = x1;
%dv1du1 = [1;0;0];
%dv1dv1 = [0;1;0];
v2 = R*x2;
%dv2du2 = R(:,1);
%dv2dv2 = R(:,2);
%dv2da1 = dR1*x2;
%dv2da2 = dR2*x2;
%dv2da3 = dR3*x2;
% distance formulas based on Schneider pp 409-412
a = sum(v1.*v1);
%dadu1 = sum(2*dv1du1.*v1);
%dadv1 = sum(2*dv1dv1.*v1);
%dadu1 = sum(2*[1;0;0].*v1);
%dadv1 = sum(2*[0;1;0].*v1);
dadu1 = 2*v1(1,:);
dadv1 = 2*v1(2,:);
dadu2 = zeros(1,n);
dadv2 = zeros(1,n);
dadt1 = zeros(1,n);
dadt2 = zeros(1,n);
dadt3 = zeros(1,n);
dada1 = zeros(1,n);
dada2 = zeros(1,n);
dada3 = zeros(1,n);
b = sum(v1.*v2);
%dbdu1 = sum(dv1du1.*v2);
%dbdv1 = sum(dv1dv1.*v2);
%dbdu1 = sum([1;0;0].*v2);
%dbdv1 = sum([0;1;0].*v2);
%dbdu2 = sum(v1.*dv2du2);
%dbdv2 = sum(v1.*dv2dv2);
dbdu1 = v2(1,:);
dbdv1 = v2(2,:);
dbdu2 = sum(v1.*R(:,1));
dbdv2 = sum(v1.*R(:,2));
dbdt1 = zeros(1,n);
dbdt2 = zeros(1,n);
dbdt3 = zeros(1,n);
%dbda1 = sum(v1.*dv2da1);
%dbda2 = sum(v1.*dv2da2);
%dbda3 = sum(v1.*dv2da3);
dbda1 = sum(v1.*(dR1*x2));
dbda2 = sum(v1.*(dR2*x2));
dbda3 = sum(v1.*(dR3*x2));
c = sum(v2.*v2);
dcdu1 = zeros(1,n);
dcdv1 = zeros(1,n);
%dcdu2 = sum(2*dv2du2.*v2);
%dcdv2 = sum(2*dv2dv2.*v2);
dcdu2 = sum(2*R(:,1).*v2);
dcdv2 = sum(2*R(:,2).*v2);
dcdt1 = zeros(1,n);
dcdt2 = zeros(1,n);
dcdt3 = zeros(1,n);
%dcda1 = sum(2*dv2da1.*v2);
%dcda2 = sum(2*dv2da2.*v2);
%dcda3 = sum(2*dv2da3.*v2);
dcda1 = sum(2*(dR1*x2).*v2);
dcda2 = sum(2*(dR2*x2).*v2);
dcda3 = sum(2*(dR3*x2).*v2);
d = xs(1:3)'*v1;
%dddu1 = xs(1:3)'*dv1du1;
%dddv1 = xs(1:3)'*dv1dv1;
dddu1 = xs(1);
dddv1 = xs(2);
dddu2 = zeros(1,n);
dddv2 = zeros(1,n);
dddt1 = v1(1,:);
dddt2 = v1(2,:);
dddt3 = v1(3,:);
ddda1 = zeros(1,n);
ddda2 = zeros(1,n);
ddda3 = zeros(1,n);
e = xs(1:3)'*v2;
dedu1 = zeros(1,n);
dedv1 = zeros(1,n);
%dedu2 = xs(1:3)'*dv2du2;
%dedv2 = xs(1:3)'*dv2dv2;
dedu2 = xs(1:3)'*R(:,1);
dedv2 = xs(1:3)'*R(:,2);
dedt1 = v2(1,:);
dedt2 = v2(2,:);
dedt3 = v2(3,:);
%deda1 = xs(1:3)'*dv2da1;
%deda2 = xs(1:3)'*dv2da2;
%deda3 = xs(1:3)'*dv2da3;
deda1 = xs(1:3)'*dR1*x2;
deda2 = xs(1:3)'*dR2*x2;
deda3 = xs(1:3)'*dR3*x2;
%f = xs(1:3)'*xr;
denom = a.*c - b.*b;
denom(denom < eps) = 1; % accounts for parallel lines
ddenomdu1 = dadu1.*c + a.*dcdu1 - 2*b.*dbdu1;
ddenomdv1 = dadv1.*c + a.*dcdv1 - 2*b.*dbdv1;
ddenomdu2 = dadu2.*c + a.*dcdu2 - 2*b.*dbdu2;
ddenomdv2 = dadv2.*c + a.*dcdv2 - 2*b.*dbdv2;
ddenomdt1 = dadt1.*c + a.*dcdt1 - 2*b.*dbdt1;
ddenomdt2 = dadt2.*c + a.*dcdt2 - 2*b.*dbdt2;
ddenomdt3 = dadt3.*c + a.*dcdt3 - 2*b.*dbdt3;
ddenomda1 = dada1.*c + a.*dcda1 - 2*b.*dbda1;
ddenomda2 = dada2.*c + a.*dcda2 - 2*b.*dbda2;
ddenomda3 = dada3.*c + a.*dcda3 - 2*b.*dbda3;
num = c.*d - b.*e;
dnumdu1 = (dcdu1.*d + c.*dddu1) - (dbdu1.*e + b.*dedu1);
dnumdv1 = (dcdv1.*d + c.*dddv1) - (dbdv1.*e + b.*dedv1);
dnumdu2 = (dcdu2.*d + c.*dddu2) - (dbdu2.*e + b.*dedu2);
dnumdv2 = (dcdv2.*d + c.*dddv2) - (dbdv2.*e + b.*dedv2);
dnumdt1 = (dcdt1.*d + c.*dddt1) - (dbdt1.*e + b.*dedt1);
dnumdt2 = (dcdt2.*d + c.*dddt2) - (dbdt2.*e + b.*dedt2);
dnumdt3 = (dcdt3.*d + c.*dddt3) - (dbdt3.*e + b.*dedt3);
dnumda1 = (dcda1.*d + c.*ddda1) - (dbda1.*e + b.*deda1);
dnumda2 = (dcda2.*d + c.*ddda2) - (dbda2.*e + b.*deda2);
dnumda3 = (dcda3.*d + c.*ddda3) - (dbda3.*e + b.*deda3);
s = num./denom;
dsdu1 = (dnumdu1.*denom - num.*ddenomdu1)./(denom.*denom);
dsdv1 = (dnumdv1.*denom - num.*ddenomdv1)./(denom.*denom);
dsdu2 = (dnumdu2.*denom - num.*ddenomdu2)./(denom.*denom);
dsdv2 = (dnumdv2.*denom - num.*ddenomdv2)./(denom.*denom);
dsdt1 = (dnumdt1.*denom - num.*ddenomdt1)./(denom.*denom);
dsdt2 = (dnumdt2.*denom - num.*ddenomdt2)./(denom.*denom);
dsdt3 = (dnumdt3.*denom - num.*ddenomdt3)./(denom.*denom);
dsda1 = (dnumda1.*denom - num.*ddenomda1)./(denom.*denom);
dsda2 = (dnumda2.*denom - num.*ddenomda2)./(denom.*denom);
dsda3 = (dnumda3.*denom - num.*ddenomda3)./(denom.*denom);
%t = (b.*d - a.*e) ./ denom;
% shortest distance between the two lines
%D = f + s.*(a.*s - b.*t - 2.*d) + t.*(c.*t - b.*s + 2.*e);
% compute lines length
r1 = sqrt(a).*s;
dr1du1 = (dadu1./(2*sqrt(a))).*s + sqrt(a).*dsdu1;
dr1dv1 = (dadv1./(2*sqrt(a))).*s + sqrt(a).*dsdv1;
dr1du2 = (dadu2./(2*sqrt(a))).*s + sqrt(a).*dsdu2;
dr1dv2 = (dadv2./(2*sqrt(a))).*s + sqrt(a).*dsdv2;
dr1dt1 = (dadt1./(2*sqrt(a))).*s + sqrt(a).*dsdt1;
dr1dt2 = (dadt2./(2*sqrt(a))).*s + sqrt(a).*dsdt2;
dr1dt3 = (dadt3./(2*sqrt(a))).*s + sqrt(a).*dsdt3;
dr1da1 = (dada1./(2*sqrt(a))).*s + sqrt(a).*dsda1;
dr1da2 = (dada2./(2*sqrt(a))).*s + sqrt(a).*dsda2;
dr1da3 = (dada3./(2*sqrt(a))).*s + sqrt(a).*dsda3;
%r2 = sqrt(b).*t;
% negative distance
%vis = r1 > 0;

% inverse depth #1 (along the rays)
% rho = 1./r1;
drhodu1 = -dr1du1./(r1.*r1);
drhodv1 = -dr1dv1./(r1.*r1);
drhodu2 = -dr1du2./(r1.*r1);
drhodv2 = -dr1dv2./(r1.*r1);
drhodt1 = -dr1dt1./(r1.*r1);
drhodt2 = -dr1dt2./(r1.*r1);
drhodt3 = -dr1dt3./(r1.*r1);
drhoda1 = -dr1da1./(r1.*r1);
drhoda2 = -dr1da2./(r1.*r1);
drhoda3 = -dr1da3./(r1.*r1);

% % inverse depth #2 (vertical to image plane)
% rim = sqrt(sum(x1(1:2,:).^2)+1);
% drimdu1 = x1(1,:)./rim;
% drimdv1 = x1(2,:)./rim;
% drimdu2 = zeros(1,n);
% drimdv2 = zeros(1,n);
% drimdt1 = zeros(1,n);
% drimdt2 = zeros(1,n);
% drimdt3 = zeros(1,n);
% drimda1 = zeros(1,n);
% drimda2 = zeros(1,n);
% drimda3 = zeros(1,n);
% %rho = rim./r1;
% drhodu1 = (drimdu1.*r1 - rim.*dr1du1) ./(r1.*r1);
% drhodv1 = (drimdv1.*r1 - rim.*dr1dv1) ./(r1.*r1);
% drhodu2 = (drimdu2.*r1 - rim.*dr1du2) ./(r1.*r1);
% drhodv2 = (drimdv2.*r1 - rim.*dr1dv2) ./(r1.*r1);
% drhodt1 = (drimdt1.*r1 - rim.*dr1dt1) ./(r1.*r1);
% drhodt2 = (drimdt2.*r1 - rim.*dr1dt2) ./(r1.*r1);
% drhodt3 = (drimdt3.*r1 - rim.*dr1dt3) ./(r1.*r1);
% drhoda1 = (drimda1.*r1 - rim.*dr1da1) ./(r1.*r1);
% drhoda2 = (drimda2.*r1 - rim.*dr1da2) ./(r1.*r1);
% drhoda3 = (drimda3.*r1 - rim.*dr1da3) ./(r1.*r1);

% Jacobian matrix

H = [drhodu1, drhodv1, drhodu2, drhodv2, ...
    drhodt1, drhodt2, drhodt3, drhoda1, drhoda2, drhoda3];   % invrese depth
%  
%  H = [dr1du1, dr1dv1, dr1du2, dr1dv2, ...
%      dr1dt1, dr1dt2, dr1dt3, dr1da1, dr1da2, dr1da3];  % depth 
 
hi = H(:);