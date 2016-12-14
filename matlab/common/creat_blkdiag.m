function b=creat_blkdiag(a,n)
% a: a square matrix
% n: dimension of output square matrix with a in the diagonal blks
% b: output matrix with a and diagonal blks
m=size(a,1);
if size(a,2)~=m;error('only for square matrices right now');end;
x=1:n;
k=m-1;
cols=[];rows=[];
for s=1:m
    i=m*x-k;
    rows=[rows repmat(i,1,m)];
    cols=[cols i];
    k=k-1;
end
cols=repmat(cols,1,m);
data=ones(size(rows));
b=full(sparse(rows,cols,data,n*m,n*m)).*repmat(a,n,n);