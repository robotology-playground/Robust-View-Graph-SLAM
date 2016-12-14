function x = hom(x)
if isempty(x)
  return
end
x(end+1,:) = 1;
return
