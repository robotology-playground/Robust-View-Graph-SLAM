function P = canonical_covariance(Y, idx)
%function P = canonical_covariance(Y, idx)
%
% INPUTS:
%   Y - information matrix
%   idx - index of states for which to recover joint covariance matrix
%
% OUTPUT:
%   P - covariance matrix for states x(idx)
%
% Compute covariance matrix for a portion of the state, x(idx). Method is
% columnwise solve, so it only incurs cost for columns P(:,idx), which may
% be much cheaper than computing all P. It may even be competitive with the
% Takahashi sparse inverse for some cases.
%
% Tim Bailey 2011.
if issparse(Y)
    I = speye(size(Y));
else
    I = eye(size(Y));
end

if nargin == 2
    P = full(Y\I(:, idx));
    P = P(idx,:);
else
    P = full(Y\I);
end
