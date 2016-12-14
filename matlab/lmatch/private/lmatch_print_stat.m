% Prints out number of tentative matches in li

function lmatch_print_stat(li)

K = size(li,1);
sli = sum(li>0,1);
n = [];
for k = 2:K
  n(k) = nnz(sli==k);
end

fprintf('(');
for k = find(n>0)
  fprintf('%i:%i ',k,n(k));
end
fprintf('any:%i',sum(n));
fprintf(')');

return
