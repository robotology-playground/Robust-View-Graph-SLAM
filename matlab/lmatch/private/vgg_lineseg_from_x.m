%vgg_lineseg_from_x   Fit line segment to a set of 2D points.
%   [u,v] = vgg_lineseg_from_x(x) fits a straight line to set of points x minimizing
%   squared orthog. distances, and finds end points of the line segment. 
%   x is double(2,?), u,v are 2-vectors.
%
%   [u,v,C] = vgg_lineseg_from_x(x) returns also covariance 3-by-3 matrix.
%   C(3,3) is number of edgels, C(:,3) is homog. coords of centroid, C(1:2,1:2) is
%   scatter matrix.
%
%   See also vgg_linesegs_from_edgestrips, vgg_xcv_segment, vgg_fit_hplane_to_x.

function [u,v,C] = vgg_lineseg_from_x(e)

C = e;
C(3,:) = 1;
C = C*C';

% l := line fitted to edgels
[U,S,V] = svd(C(1:2,1:2)-C(1:2,3)*C(1:2,3)'/C(3,3),0);
a = U(:,end)';
l = [a -a*C(1:2,3)/C(3,3)];

% end points
d = vgg_contreps(l(1:2))'*e;
ab = [-l(2) l(1)];
cl = vgg_contreps(l);
[dummy,imin] = min(d);
[dummy,imax] = max(d);
u = cl*[ab -ab*e(:,imin)]';  u = u(1:2)/u(3);
v = cl*[ab -ab*e(:,imax)]';  v = v(1:2)/v(3);

return