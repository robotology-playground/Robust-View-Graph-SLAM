function s = cellsizes(C)
for k = 1:prod(size(C))
  s(:,k) = size(C{k})';
end
return