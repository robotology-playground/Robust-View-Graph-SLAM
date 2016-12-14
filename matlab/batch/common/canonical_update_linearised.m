function [y, Y, w, v, S] = canonical_update_linearised(y, Y, model, norm, z, R, xs, Hs, idx, islog, varargin)
%[y, Y, w, v, S] = canonical_update_linearised(y, Y, model, norm, z, R, xs, Hs, idx, islog, varargin)
%
% INPUTS:
%   y,Y - information form Gaussian
%   model - function handle (eg., @h or 'h') for observation model
%   norm - function handle for innovation normalisation v = norm(z-zpred)
%   z,R - observation and its covariance
%   xs - linearisation point for model
%   Hs - Jacobian for linearised observation model
%   idx - index such that z = h(x(idx), ...) 
%   islog - flag for update weight (w), see below
%   ... - extra variables in z = h(x(idx), ...)
%
% OUTPUTS:
%   y,Y - updated Gaussian
%   w - weight of update p(z) (ie., Bayes evidence)
%   v,S - innovation and its covariance
%
% Tim Bailey, 2011.

% EXPENSIVE CHECK: Hs*xs below may cause a problem if there is a
% discontinuity between x_mean and xs, since x_mean is encoded in y.
if debug_conditional(3) && rank(full(Y)) == size(Y,1)
    xm = Y\y;
    P = canonical_covariance(Y, idx);
    check_for_discontinuity(xm(idx)-xs, P);
end

% Calculate weight of update
% Caution: expensive. We have to recover mean and covariance of idx-states
if nargout > 2
    [v, S] = canonical_innovation_linearised(y, Y, model, norm, z, R, xs, Hs, idx, varargin{:});
    w = moment_evaluate(v, S, 0, islog);
end

% Update information matrix
HtRi = Hs'/R;
Y(idx,idx) = Y(idx,idx) + force_symmetry(HtRi*Hs);

% Update information vector
zs = feval(model, xs, varargin{:});
v = z - zs;
% %%% adde by Tariq, using rotation matrix to compute angular difference
% Rzst = w2R(zs(4:6))';
% Rz = w2R(z(4:6));
% Rv = Rzst*Rz;
% v(4:6) = R2w(Rv)'; 
% %%%

if ~isempty(norm)
    v = feval(norm, v); 
    % FIXME: Is this the correct way to address discontinuous models here,
    % since z-zs is *not* the innovation (= z - zpred)? [Note: zpred = zs
    % + Hs(x_mean - xs)]. Is there a counter-example? I think there may be
    % a problem if x_mean - xs has a discontinuity.  
end
y(idx) = y(idx) + HtRi*(v + Hs*xs);