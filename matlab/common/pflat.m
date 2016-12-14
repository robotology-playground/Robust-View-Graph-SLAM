function y = pflat(x)

%y = x./repmat(x(end,:),[size(x,1) 1]);
y = x./(ones(size(x,1),1)*x(end,:));