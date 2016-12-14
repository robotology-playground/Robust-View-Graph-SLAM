function e = ellipse_sigma(x, P, nsig, N, type)
%function e = ellipse_sigma(x, P, nsig, N, type)
%
% INPUTS:
%   x, P - mean and covariance of a 2-D Gaussian
%   nsig - number of sigmas for ellipsoid bound
%   N - number of lines in polyline (default 60)
%
% OUTPUT:
%   e = points of ellipse polyline
%
% Tim Bailey 2006, modified 2011.

if nargin>3, inc = 2*pi/N; else inc = pi/30; end
if nargin<5, type = 1; end

switch type
    case 1 
        r = chol(P)'; % note, must be lower triangular
    case 2
        r = sqrtm(P);
    case 3
        [v,d] = eig(P);
        r = v*sqrt(d);        
end

phi = [0:inc:2*pi 0];
a = nsig * r * [cos(phi); sin(phi)];

e(1,:) = a(1,:) + x(1);
e(2,:) = a(2,:) + x(2);
