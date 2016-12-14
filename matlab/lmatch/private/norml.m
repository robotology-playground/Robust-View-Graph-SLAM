% l = norml(l)  Multiplies hyperplane l by scalar so that for each n, norm(l(1:end-1,n))==1. 
function l = norml(l)
l = l./(sqrt(sum(l(:,1:end-1).^2,2))*ones(1,size(l,2)));
return