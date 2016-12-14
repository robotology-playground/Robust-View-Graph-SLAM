function y = normalise_points(x, K, d)

% Computes the normalization of projective points x so the last coordinate
% becomes 1
% Input:
% x     is N points in image coordinates (2xN)
% K     is the camera matrix containing intrinsic parameters
% Output:
% y     is N points in homogenous and normalised coordinates (3xN) 
%       with last row ones
%
%   Tariq Abuhashim

if nargin<3;
    d = [];
end

[m, n] = size(x);
rot = 0;

% check if none of the dimensions is 2 or 3
if (m > 3) || (m < 2 && n > m)
    error(' points should be in column vectors 2xN or 3xN');
end

% check the orientation of the input vector and make sure each element is a
% column vector (ie, size(x, 1)==2 or size(x, 1)==3)
%if ((m > n) && (n > 3)) || ((n > m) && (m < 3)); x = x'; rot = 1; end

% extend if needed
if size(x,1) < 3; x = pextend(x); end

% calibrate using the intrinsic camera matrix(K)
y = K \ x;

% flatten, so that 
y = pflat(y);

% undistort points too
if ~isempty(d)
    y = cv.undistortPoints(x, K, d);
end

% fix the orientation of the output dimensions to match those of the input
%if rot; y = y(1:2, :)'; end