% L = lineseg_merge(K [,x])  Merges several line segments K, optionally going through point x.
%
% K ... lineseg (1,N)
% x ... double (3,1)
% L ... lineseg

function L = lineseg_merge(K,x)

if length(K)==1
  L = K;
  return
end

K = K(:)';
M = length(K);
L = K(1);

% new covariance matrix
L.s = sum([K.s],2);
s = vgg_vech(L.s);

% l := line fitted to edgels
if nargin == 1 % unconstrained fit
  [U,S,V] = svd(s(1:2,1:2)-s(1:2,3)*s(1:2,3)'/s(3,3));
  a = [0 1]*U';
  l = norml([a -a*nhom(s(:,3))]);
else % l is constrained to go through x
  l = norml(s(3,:)*vgg_contreps(x));
end

% find new end points
e = [[K.u] [K.v]];
d = vgg_contreps(l(1:2))'*e;
ab = [-l(2) l(1)];
cl = vgg_contreps(l);
L.u = nhom(cl*[ab -ab*e(:,argmin(d))]');
L.v = nhom(cl*[ab -ab*e(:,argmax(d))]');

% orient resulting segment to accordance with the first source segment
if (L.v-L.u)'*(K(1).v-K(1).u) < 0
  aux = L.v;
  L.v = L.u;
  L.u = aux;
end

return