function [dim, flag] = test_point_feature_dim(p)
%
% [dim, flag] = test_point_feature_dim(p)
% This function tests the properties of the input image point feature
% and outputs a flag as follows:
% 
% dim:
% 0 - if a row vector
% 1 - if a column vector
% flag:
% 0 - if not in homogeneous coordinates
% 1 - if in homogeneous coordinates
% 2 - if in homogeneous coordinates and is normalised
%
% Tariq Abuhashim - 2015
%

% column vector test
if size(p, 1) == 2 && size(p, 2);
if size(p, 1) > 1 && size(p, 1) < 4;
    if size(p, 1) == 3 && sum(p(3,:))==size(p,2); % homogeneous coordinates test
        flag = 1;
    end
end


