function x = nhom(x)
if isempty(x)
  x = []; 
  return; 
end
d = size(x,1) - 1;
x = x(1:d,:)./(ones(d,1)*x(end,:));
return
