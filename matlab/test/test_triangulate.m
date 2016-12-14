function [xf, r1, vis] = test_triangulate(Rr, xr, x1, x2)

% function [xf, J, D, r1, r2, vis] = test_triangulate(Rr, xr, x1, x2)
% computes triangulation and its parameters
% Rr : rotation
% xr : translation
% x1 : feature coordinates in the first image
% x2 : feature coordinates in the second image
%
% to test the analytical Jacobian against numerical Jacobian:
% Rr = eye(3); xr = rand(3,1); x1 = rand(2,1); x2 = x1 + [.1;0];
% [xf, J] = test_triangulate(Rr, xr, x1, x2);
% J1 = numerical_jacobian_i(@test_triangulate, [], 3, [], Rr, t, x1, x2);
% J2 = numerical_jacobian_i(@test_triangulate, [], 4, [], Rr, t, x1, x2);
% error = full(J) - [J1 J2]
%
% Tariq Abuhashim - 2015
%

if size(x1,1) < 3; x1 = pextend(x1); end;
if size(x2,1) < 3; x2 = pextend(x2); end;

% calculate lines
v1 = x1;
%v1 = v1./repmat(sqrt(sum(v1.^2)), [3, 1]);
v2 = Rr*x2;
%v2 = v2./repmat(sqrt(sum(v2.^2)), [3, 1]);

% distance formulas based on Schneider pp 409-412
a = sum(v1.*v1);
b = sum(v1.*v2);
c = sum(v2.*v2);
d = xr'*v1;
e = xr'*v2;
%f = xr'*xr;
denom = a.*c - b.*b;
denom(denom < eps) = 1; % accounts for parallel lines
num = c.*d - b.*e;
s = num./denom;
%t = (b.*d-a.*e)./denom;

% shortest distance between the two lines
%D = f + s.*(a.*s-b.*t-2.*d)+t.*(c.*t-b.*s+2.*e);

% compute lines length
r1 = sqrt(a).*s;
%r2 = sqrt(b).*t;

% negative distance
vis = r1>0;

% add uncertainty model:
% x = (u*sigma(d)^2 + d*sigma(u)^2) / (sigma(d)^2*sigma(u)^2)

% 3d points in P1 coordinates
rim = sqrt(sum(x1(1:2, :).^2) + 1);
z = r1./rim;
xf = zeros(3,size(x1,2));
xf(1,:) = x1(1,:).*z;
xf(2,:) = x1(2,:).*z;
xf(3,:) = x1(3,:).*z;

% %%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Perpendicular to the lines v1 cross v2
% n = [v1(2,:).*v2(3,:)-v1(3,:).*v2(2,:); ...
%      v1(3,:).*v2(1,:)-v1(1,:).*v2(3,:); ...
%      v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:)];
% n = n./repmat(sqrt(sum(n.^2)),[3 1]);
%
% % Shortest vector between the lines
% d = repmat(-xr'*n,[3 1]).*n; % using repmat
%
% % Point on l2 closest to l1
% % Plane containing n and l1
% npi = [v1(2,:).*n(3,:)-v1(3,:).*n(2,:); ...
%        v1(3,:).*n(1,:)-v1(1,:).*n(3,:); ...
%        v1(1,:).*n(2,:)-v1(2,:).*n(1,:)];
%
% % Cutting with l2
% lambda = -xr'*npi./(sum(npi.*v2));
% p2 = repmat(xr,[1 size(v2,2)])+repmat(lambda,[3 1]).*v2;  % using repmat
%
% % Mid-point
% xf = p2+d/2;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% P1 = eye(3,4); P2 = [Rr' Rr'*xr];
% for i=1:size(x1,2)
%     %Select point
%     sx1 = x1(:,i);
%     sx2 = x2(:,i);
%     A1 = sx1(1,1).*P1(3,:) - P1(1,:);
%     A2 = sx1(2,1).*P1(3,:) - P1(2,:);
%     A3 = sx2(1,1).*P2(3,:) - P2(1,:);
%     A4 = sx2(2,1).*P2(3,:) - P2(2,:);
%     %Set A matric, and find SDV(A)
%     A = [A1;A2;A3;A4];
%     [U,D,V] = svd(A);
%     %Point in 3D space is the last column of V
%     X_temp = V(:,4);
%     X_temp = X_temp ./ repmat(X_temp(4,1),4,1);
%     xf(:,i) = X_temp;
% end
% r1 = [];
% r2 = [];
% D = [];
% vis = ones(1, size(xf,2));
%
% %%%%%%%%%%%%%%%%%%%%%%%%%