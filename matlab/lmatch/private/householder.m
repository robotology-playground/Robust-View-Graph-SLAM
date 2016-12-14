% R = householder(x,y)  Given x, y UNIT vectors, R is orthogonal such that R*x=y, R==R'.

function R = householder(x,y)

% choose corect sign of y to reduce roundoff error
if x'*y < 0
  u = x-y;
  R = eye(size(u,1)) - 2*u*u'/(u'*u);
else
  u = x+y;
  R = -eye(size(u,1)) + 2*u*u'/(u'*u);
end

return