function C = apply_nominal_R(C, R)
%function C = apply_nominal_R(C, R)
%
% INPUTS:
%   C - set of constraints C(i,j).z between poses i and j. If there is no
%   constraint between two poses, C(i,j).z is empty.
%   R - nominal constraint covariance, to be applied to all constraints
%
% OUTPUT:
%   C - constraints now including nominal R
%
% Tariq Abuhashim, 2015.

if nargin == 1
    config_rswitch_v2;
    R = diag([SIGMA_t,SIGMA_t,SIGMA_t,SIGMA_a,SIGMA_a,SIGMA_a].^2);
    %R = inv(C(i).Y);
end

for i = 1:length(C)
    C(i).R = R;
end