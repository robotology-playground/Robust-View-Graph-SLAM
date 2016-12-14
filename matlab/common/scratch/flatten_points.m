function y = flatten_points(x)

% Computes the normalization of projective points x so the last coordinate
% becomes 1
% Input:
% x     is N points in image coordinates (3xN)
% Output:
% y     is N points in homogenous and normalised coordinates (3xN) 
%       with last row ones
%
%   Tariq Abuhashim



[rows,npts] = size(x);
y = x;

% Find the indices of the points that are not at infinity
finiteind = find(abs(x(rows,:)) > eps);
if length(finiteind) ~= npts
	warning('Some points are at infinity');
end

% Normalise points not at infinity
x = x(:, finiteind);
alpha = x(rows, :);
a = size(x, 1);
y(:, finiteind) = x./(ones(a, 1)*alpha);

% the old version
%a = size(x,1);
%alpha = x(a,:);
%y = x./(ones(a,1)*alpha);