function Q = transform_posdef(P, F)
%function Q = transform_posdef(P, F)
%
% Perform transform Q = F*P*F', where P is a symmetric positive definite
% matrix. Guarantees Q is symmetric positive definite. 
%
% Caveat: This function is expensive if P is large. In such cases, one may
% prefer merely to compute Q naively as F*P*F' and then enforce symmetry,
% Q = (Q+Q')/2. 
%
% Tim Bailey 2011.

R = chol(P);    
RF = R*F';
Q = RF'*RF; % nice, guarantees symmetry and pos-def
