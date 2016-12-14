% x = normx(x)  Normalize MxN matrix so that norm of each its column is 1.
function x = normx(x)
if ~isempty(x)
  x = x./(ones(size(x,1),1)*sqrt(sum(x.*x)));
end
return