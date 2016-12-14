function y = extend_points(x)

% Computes the normalization of projective points x so the last coordinate
% becomes 1
% Input:
% x     is N points in non-homogenous coordinates (2xN)
% Output:
% y     is N points in homogenous normalised coordinates (3xN)-last row 1s
%
%   Tariq Abuhashim



a = size(x,2);
y=[x; ones(1,a)];

