% [Dnear,k0] = extreme_view(K,k,sgn,D)  Finds nearest or farthest view of view set K to view k.
%
% K ... vector of view indices, set of views
% sgn ... +1 for nearest view, -1 for farthest view
% D ... double(#view,#views), distance table of pairs of views
% Dnear ... distance of the nearest view
% k0 ... its index (element of K)

function [Dnear,k0] = extreme_view(K,k,sgn,D)

Dnear = +Inf*sgn;
for l = K(:)'
  d = D(l,k);
  if d*sgn < Dnear*sgn
    k0 = l;
    Dnear = d;
  end
end

return