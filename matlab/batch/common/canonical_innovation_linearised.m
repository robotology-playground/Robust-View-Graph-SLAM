function [v,S] = canonical_innovation_linearised(y, Y, model, norm, z, R, xs, Hs, idx, varargin)
%function [v,S] = canonical_innovation_linearised(y, Y, model, norm, z, R, xs, Hs, idx, ...)
%
% INPUTS:
%   y,Y - information form Gaussian
%   model - function handle (eg., @h or 'h') for observation model
%   norm - function handle for innovation normalisation v = norm(z-zpred)
%   z,R - observation and its covariance
%   xs - linearisation point for model
%   Hs - Jacobian for linearised observation model
%   idx - index such that z = h(x(idx), ...) 
%   ... - extra variables in z = h(x(idx), ...)
%
% OUTPUTS:
%   v,S - innovation and its covariance
%
% Tim Bailey 2011.

% Recover moments
x = Y \ y;
P = canonical_covariance(Y, idx);
if debug_conditional(1)
    check_for_discontinuity(x(idx) - xs, P);
end

% Compute innovation covariance
S = transform_posdef(P, Hs) + R;

% Compute innovation
zpred = feval(model, xs, varargin{:}) + Hs * (x(idx) - xs); % predicted observation
v = z - zpred;
if ~isempty(norm)
    v = feval(norm, v);
end
