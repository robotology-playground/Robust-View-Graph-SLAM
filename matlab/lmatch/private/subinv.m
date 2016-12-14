% subinv  See subtx.

function P = subinv(P)

p = P(1:end-1,1:end-1)';
P = [p -p*P(1:end-1,end)/P(end,end); zeros(1,size(P,1)-1) 1];

return