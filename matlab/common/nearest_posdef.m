function out = nearest_posdef(in)

% Author: MB (marco.buchmann@ecb.int)
% Date: Dec 2011
%
% The function computes a positive semi-definite matrix 'out' that is 
% closest to a symmetric covariance matrix 'in' which is not positive
% semi-definite. 

% Check/process input
if sum(sum(in-in'))~=0
    error('The input matrix must be symmetric')
end
if sum(diag(in)>0)~=length(diag(in))
    error('The input matrix shall have positive elements on its diagonal')
end
if min(eig(in))>=0
    out = in;
    disp('The covariance matrix is already positive semi-definite')
    return
end
global Size; Size = size(in);
% Process upper triangular part of input and make it global
f = reshape(triu(in)',[],1); f(f==0) = []; global F; F = f;
% Set zero lower bounds for variances
lbnd = NaN(length(f),1); diagin = diag(in);
for i=1:length(diagin)
    for j=1:length(f)
        if f(j)==diagin(i)
            lbnd(j) = 0;
        elseif f(j)~=diagin(i) && lbnd(j)~=0
            lbnd(j) = -Inf;
        end
    end
end
% Invoke nonlinear optimization
back = fmincon(@distF,f,[],[],[],[],lbnd,[],@ineq,optimset('Algorithm','active-set','MaxFunEvals',20000,'MaxIter',20000,'Display','off','PlotFcns',@optimplotfirstorderopt));
close
% Bring output back to symmetric matrix format
a = zeros(size(in)); ct = 1; rc = 1; cc = 1;
for i=1:length(back)
    a(rc,cc) = back(i); cc = cc+1;
    if cc>Size(2)
        cc = ct+1; rc = rc+1; ct = ct+1;
    end
end
out_pre = a+a';
for i=1:Size(1)
    out_pre(i,i) = a(i,i);
end
% Should small negative eigenvalues remain, correct them
[ve,va] = eig(out_pre); va(va<0) = 0; out = ve*va*ve';
% Visualize difference between in and out
diffF = in-out; imagesc(diffF(1:size(diffF,1),1:size(diffF,2))); colorbar; colormap(copper)
set(gcf,'Color',[.97 .97 .97]); title('Difference between matrix input and output','FontSize',12)
end

function out = distF(in) % distance
global F;
out = sum((F-in).^2);
end

function [c,ceq] = ineq(in) % inequality constraints (eigenvalues larger zero)
global Size;
a = zeros(Size); ct = 1; rc = 1; cc = 1;
for i=1:length(in)
    a(rc,cc) = in(i); cc = cc+1;
    if cc>Size(2)
        cc = ct+1; rc = rc+1; ct = ct+1;
    end
end
b = a+a';
for i=1:Size(1)
    b(i,i) = a(i,i);
end
[ve,va] = eig(b); c = -diag(va); ceq = [];
end
