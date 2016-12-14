function [rho, J, vis] = test_triangulate_inverse_depth(p1, p2, xs, varargin)
%[rho, r1, vis] = test_triangulate_inverse_depth(p1, p2, xs)
%
% computes inverse depth triangulation and it jacobian
% p1 : feature coordinates in the first image
% p2 : feature coordinates in the second image
% xs : camera motion parameters {t,r}
%
% to test the analytical Jacobian against numerical Jacobian:
% x1 = rand(2,1); x2 = rand(2,1); xs = rand(6,1);
% [rho, J] = test_triangulate_inverse_depth(x1, x2, xs);
% J1 = numerical_jacobian_i(@test_triangulate_inverse_depth, [], 1, [], x1, x2, xs);
% J2 = numerical_jacobian_i(@test_triangulate_inverse_depth, [], 2, [], x1, x2, xs);
% J3 = numerical_jacobian_i(@test_triangulate_inverse_depth, [], 3, [], x1, x2, xs);
% error = full(J) - [J1 J2 J3]
%
% Tariq Abuhashim, 2015.
%
% iCub - Koroibot

if size(p1,1) < 3; p1 = pextend(p1); end;
if size(p2,1) < 3; p2 = pextend(p2); end;

% calculate lines
v1 = p1;
%v1 = v1./repmat(sqrt(sum(v1.^2)), [3, 1]);
v2 = w2R(xs(4:6))*p2;
%v2 = v2./repmat(sqrt(sum(v2.^2)), [3, 1]);

% distance formulas based on Schneider pp 409-412
a = sum(v1.*v1);
b = sum(v1.*v2);
c = sum(v2.*v2);
d = xs(1:3)'*v1;
e = xs(1:3)'*v2;
%f = xr'*xr;
denom = a.*c - b.*b;
denom(denom < eps) = 1; % accounts for parallel lines
num = c.*d - b.*e;
s = num./denom;

% compute lines length
r = sqrt(a).*s;

% negative distance
vis = r > 0;

% inverse depth #1 (along the rays)
rho = 1./r;  % inverse depth
%rho = r;  % depth 
% the depth would match test_triangulate_jacobian_inverse_depth output
% however, for estimation purposes, test_triangulate_inverse_depth outputs
% the inverse depth value, but test_triangulate_jacobian_inverse_depth
% outputs the jacobian of the depth, instead.

% % inverse depth #2 (vertical to image plane)
% rim = sqrt(sum(p1(1:2,:).^2)+1);
% rho = rim./r;

if nargout > 1;
    J = test_triangulate_jacobian_inverse_depth(p1, p2, xs, varargin{:});
end
