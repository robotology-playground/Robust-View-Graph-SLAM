function P = force_symmetry(P)
% function P = force_symmetry(P)
%
% INPUT:
%   P - square matrix, possibly not symmetric
%
% OUTPUT:
%   P - symmetric matrix
%
% Force matrix to be symmetric, P == P', by averaging off-diagonal terms.
%
% Tim Bailey 2009.

P = (P+P') / 2;
