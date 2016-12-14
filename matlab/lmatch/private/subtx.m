% subtx  Orthonormal parametrization of a subspace of a linear space.
%
% [p,x0] = subtx(S) returns transformation of restriction/projection from
% linear space V to linear subspace S. The transformation preserves metric
% in the embedded metric space.
% .. or P = subtx(S)
%
% P = subtx(S) is orthogonalization of S such that S(end,1:end)==0. P = [p x0; 0..0 a].
%
% S ... K-dimensional subspace of D-dimensional linear space (K<D).
%       The subspace is defined by :-
%    - either, if size(S)=[D K], by join of K points (points are columns of S)
%    - or, if size(S)=[K D], by meet of K hyperplanes (hyperplanes are rows of S)
% p ... double (D-1,K), affine projection matrix, p'*p==eye(K)
% x0 ... double(D-1,1), affine origin, orthogonal projection of origin to S
%
% If size(S)=[D K], the sign of P is chosen such that det(R)>0 where subtx(S)*R = S.
%
% Projection equations (x and X in V, m and M in S) :-
%  - in non-homogeneous points: x = p*m + x0, m = p'*(x - x0), 
%  - in homogeneous points: X = P*M, M = subinv(P)*X
%   (subinv(P) does just simple rearrangment of P, for convenenience)
%
% Works only for finite subspaces (not eg 3D lines at infinity).

% - Closest points in two different subspaces (eg, on two disjoint 3D lines).
%   If the subspaces are S1 and S2, and [p1,x01]=subtx(S1), [p2,x02]=subtx(S2),
%   then the closest points m1, m2 are given by 
%        [m1;m2] = [p1 -p2]\(x02-x01).


function [P,x0] = subtx(S)

[N K] = size(S);

% P := unitary matrix with columns spanning the subspace
if N>=K
  if nnz(abs(eye(K)-S'*S)>100*eps)>0 % check whether S is already not orthogonal
    [P,s,U] = svd(S,0);
  end
else
  P = null(S);
  [N K] = size(P);
end

% make P(end,1:end-1) = 0
P = P*householder(normx(P(end,:)'),[zeros(K-1,1);1]);

% signs
if size(S,1) > size(S,2)
  P = P * sign(det(subinv(P)*S)*P(end,end));
end

% make P(end,end)==1
P(:,end)=P(:,end)/P(end,end);

if nargout>1
  x0 = P(1:end-1,end);
  P = P(1:end-1,1:end-1);
end

return





% choose sign of P such that wedge{P} ==+ wedge{S}
if size(S,1)>size(S,2) % otherwise, I dont know how to check signs in computationally feasible way
  if N > 10
    warning('Sign of wedge{P} ignored - N too large.');
  else
    for i = nchoosek(1:N,K)'
      dP = det(P(i,:));
      if abs(dP) > 1e-6 % determinant is large enough
	P = P*sign(dP*det(S(i,:)));
      end
      break
    end
  end
end





return