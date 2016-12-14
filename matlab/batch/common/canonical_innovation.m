function [v,S] = canonical_innovation(y,Y, z,R,H, idx)
%function [v,S] = canonical_innovation(y,Y, z,R,H, idx)
%
% INPUTS:
%   y,Y - information form Gaussian
%   z,R - observation and its covariance
%   H - Jacobian for linear observation model
%   idx - model has form, z = H*x(idx) + r, where r ~ N(0,R)
%
% OUTPUTS:
%   v,S - innovation and its covariance
%
% Warning, this function involves recovering {x(idx), P(idx,idx)}, which is
% relatively expensive.
%
% Tim Bailey 2011.

x = Y \ y;
P = canonical_covariance(Y, idx);
v = z - H*x(idx);
S = transform_posdef(P, H) + R;

if debug_conditional(1)
    check_for_discontinuity(v, S);
end
