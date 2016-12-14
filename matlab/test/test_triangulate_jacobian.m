function H = test_triangulate_jacobian(Rr, xr, x1, x2)
% 
% function H = test_triangulate_jacobian(Rr, xr, x1, x2)
% computes Jacobian of features triangulation
% Rr : rotation
% xr : translation
% x1 : feature coordinates in the first image
% x2 : feature coordinates in the second image
%
% to test against numerical Jacobian:
% R = eye(3); t = rand(3,1); x1 = rand(2,1); x2 = x1 + [.1;0]; 
% J = test_triangulate_jacobian(R, t, x1, x2);
% J1 = numerical_jacobian_i(@test_triangulate, [], 3, [], R, t, x1, x2);
% J2 = numerical_jacobian_i(@test_triangulate, [], 4, [], R, t, x1, x2); 
% error = full(J) - [J1 J2]
%
% Tariq Abuhashim - 2015
%

% Allocate storage for triplet form
N = size(x1, 2);
Hii = zeros(12*N, 1); % 12 diagonal terms for each point
for i = 1:N
    ih = getindex_point_triplet(i);
    hi = compute_derivatives(Rr, xr, x1(:, i), x2(:, i));
    Hii(ih) = hi;
end
[ii, jj] = compute_block_diagonal_triplet_indices(N);
Hii = sparse(ii, jj, Hii, 3*N, 4*N);
H = Hii;

%
%
function c = make_col(m)
c = m(:);
%
%
function idx = getindex_point_triplet(i)
idx = (1:12) + (i-1)*12;
%
%
function [i, j] = compute_block_diagonal_triplet_indices(N)
% Compute indices of the triplet-form block-diagonal terms
% 3D points
i1 = zeros(3, N);
i1(:) = 1:3*N;
i1 = make_col(repmat(i1, 4, 1));
j1 = make_col(reprow(1:4*N, 3));
i = i1;
j = j1;
%
%
function hi = compute_derivatives(Rr, xr, x1, x2)
if size(x1, 1) < 3; x1 = pextend(x1); end
if size(x2, 1) < 3; x2 = pextend(x2); end
% calculate the terms
v1 = x1;
n = size(v1,2);
dv1du1 = repmat([1;0;0], 1, n);
dv1dv1 = repmat([0;1;0], 1, n);
v2 = Rr*x2;
dv2du2 = repmat(Rr(:,1), 1, n);
dv2dv2 = repmat(Rr(:,2), 1, n);
% distance formulas based on Schneider pp 409-412
a = sum(v1.*v1);
dadu1 = sum(dv1du1.*v1 + v1.*dv1du1);
dadv1 = sum(dv1dv1.*v1 + v1.*dv1dv1);
dadu2 = zeros(1, n);
dadv2 = zeros(1, n);
b = sum(v1.*v2);
dbdu1 = sum(dv1du1.*v2);
dbdv1 = sum(dv1dv1.*v2);
dbdu2 = sum(v1.*dv2du2);
dbdv2 = sum(v1.*dv2dv2);
c = sum(v2.*v2);
dcdu1 = zeros(1, n);
dcdv1 = zeros(1, n);
dcdu2 = sum(dv2du2.*v2 + v2.*dv2du2);
dcdv2 = sum(dv2dv2.*v2 + v2.*dv2dv2);
d = xr'*v1;
dddu1 = xr'*dv1du1;
dddv1 = xr'*dv1dv1;
dddu2 = zeros(1,n);
dddv2 = zeros(1,n);
e = xr'*v2;
dedu1 = zeros(1,n);
dedv1 = zeros(1,n);
dedu2 = xr'*dv2du2;
dedv2 = xr'*dv2dv2;
%f = xr'*xr;
denom = a.*c - b.*b;
denom(denom < eps) = 1; % accounts for parallel lines
ddenomdu1 = dadu1.*c + a.*dcdu1 - 2*b.*dbdu1;
ddenomdv1 = dadv1.*c + a.*dcdv1 - 2*b.*dbdv1;
ddenomdu2 = dadu2.*c + a.*dcdu2 - 2*b.*dbdu2;
ddenomdv2 = dadv2.*c + a.*dcdv2 - 2*b.*dbdv2;
num = c.*d - b.*e;
dnumdu1 = (dcdu1.*d + c.*dddu1) - (dbdu1.*e + b.*dedu1);
dnumdv1 = (dcdv1.*d + c.*dddv1) - (dbdv1.*e + b.*dedv1);
dnumdu2 = (dcdu2.*d + c.*dddu2) - (dbdu2.*e + b.*dedu2);
dnumdv2 = (dcdv2.*d + c.*dddv2) - (dbdv2.*e + b.*dedv2);
s = num./denom;
dsdu1 = (dnumdu1.*denom - num.*ddenomdu1)./(denom.*denom);
dsdv1 = (dnumdv1.*denom - num.*ddenomdv1)./(denom.*denom);
dsdu2 = (dnumdu2.*denom - num.*ddenomdu2)./(denom.*denom);
dsdv2 = (dnumdv2.*denom - num.*ddenomdv2)./(denom.*denom);
%t = (b.*d - a.*e) ./ denom;
% shortest distance between the two lines
%D = f + s.*(a.*s - b.*t - 2.*d) + t.*(c.*t - b.*s + 2.*e);
% compute lines length
r1 = sqrt(a).*s;
dr1du1 = ( dadu1./(2*sqrt(a)) ).*s + sqrt(a).*dsdu1;
dr1dv1 = ( dadv1./(2*sqrt(a)) ).*s + sqrt(a).*dsdv1;
dr1du2 = ( dadu2./(2*sqrt(a)) ).*s + sqrt(a).*dsdu2;
dr1dv2 = ( dadv2./(2*sqrt(a)) ).*s + sqrt(a).*dsdv2;
%r2 = sqrt(b).*t;
% negative distance
%vis = r1 > 0;
% 3d points in P1 coordinates
rim = sqrt(sum(x1(1:2, :).^2) + 1);
drimdu1 = x1(1, :)./rim;
drimdv1 = x1(2, :)./rim;
drimdu2 = zeros(1, n);
drimdv2 = zeros(1, n);
z = r1./rim;
dzdu1 = (dr1du1.*rim - r1.*drimdu1)./(rim.*rim);
dzdv1 = (dr1dv1.*rim - r1.*drimdv1)./(rim.*rim);
dzdu2 = (dr1du2.*rim - r1.*drimdu2)./(rim.*rim);
dzdv2 = (dr1dv2.*rim - r1.*drimdv2)./(rim.*rim);
%xf = x1.*repmat(z,3,1);
dXdu1 = dv1du1.*repmat(z, 3, 1) + x1.*repmat(dzdu1, 3, 1);
dXdv1 = dv1dv1.*repmat(z, 3, 1) + x1.*repmat(dzdv1, 3, 1);
dXdu2 = x1.*repmat(dzdu2, 3, 1);
dXdv2 = x1.*repmat(dzdv2, 3, 1);
% Jacobian matrix
H = [dXdu1, dXdv1, dXdu2, dXdv2]; hi = H(:);
% hi = zeros(12, 1);
% hi(1) = dXdu1(1,:);
% hi(2) = dXdu1(2,:);
% hi(3) = dXdu1(3,:);
% hi(4) = dXdv1(1,:);
% hi(5) = dXdv1(2,:);
% hi(6) = dXdv1(3,:);
% hi(7) = dXdu2(1,:);
% hi(8) = dXdu2(2,:);
% hi(9) = dXdu2(3,:);
% hi(10) = dXdv2(1,:);
% hi(11) = dXdv2(2,:);
% hi(12) = dXdv2(3,:);


% OLD CODE
% ===========

% % fill-up the jacobian matrix
% lastentry = 0;
% % dXdu1
% row(lastentry+(1:n)) = (1:3:3*n)';
% col(lastentry+(1:n)) = (1:4:4*n)';
% temp=dxfdu1(1,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dXdv1
% row(lastentry+(1:n))=(1:3:3*n)';
% col(lastentry+(1:n))=(2:4:4*n)';
% temp=dxfdv1(1,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dXdu2
% row(lastentry+(1:n))=(1:3:3*n)';
% col(lastentry+(1:n))=(3:4:4*n)';
% temp=dxfdu2(1,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dXdv2
% row(lastentry+(1:n))=(1:3:3*n)';
% col(lastentry+(1:n))=(4:4:4*n)';
% temp=dxfdv2(1,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dYdu1
% row(lastentry+(1:n))=(2:3:3*n)';
% col(lastentry+(1:n))=(1:4:4*n)';
% temp=dxfdu1(2,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dYdv1
% row(lastentry+(1:n))=(2:3:3*n)';
% col(lastentry+(1:n))=(2:4:4*n)';
% temp=dxfdv1(2,:);
% data(lastentry+(1:n))=temp';
% lastentry= lastentry+n;
% % dYdu2
% row(lastentry+(1:n))=(2:3:3*n)';
% col(lastentry+(1:n))=(3:4:4*n)';
% temp=dxfdu2(2,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dYdv2
% row(lastentry+(1:n))=(2:3:3*n)';
% col(lastentry+(1:n))=(4:4:4*n)';
% temp=dxfdv2(2,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dZdu1
% row(lastentry+(1:n))=(3:3:3*n)';
% col(lastentry+(1:n))=(1:4:4*n)';
% temp=dxfdu1(3,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dZdv1
% row(lastentry+(1:n))=(3:3:3*n)';
% col(lastentry+(1:n))=(2:4:4*n)';
% temp=dxfdv1(3,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dZdu2
% row(lastentry+(1:n))=(3:3:3*n)';
% col(lastentry+(1:n))=(3:4:4*n)';
% temp=dxfdu2(3,:);
% data(lastentry+(1:n))=temp';
% lastentry=lastentry+n;
% % dZdv2
% row(lastentry+(1:n))=(3:3:3*n)';
% col(lastentry+(1:n))=(4:4:4*n)';
% temp=dxfdv2(3,:);
% data(lastentry+(1:n))=temp';
% % fill-in the sparse jacobian
% data(isnan(data))=0;
% %data(data>1e4)=0;
% H = sparse(row, col, data, 3*n, 4*n);