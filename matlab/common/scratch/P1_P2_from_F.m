function [P1, P2] = P1_P2_from_F(F)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT USE THIS FUNCTION, BECAUSE THE SOLUTION HAS SIGN AND DIRECTION  %
% AMBIGUITIES. THERE SHOULD BE FOUR SOLUTIONS, NOT ONLY ONE GIVEN F or E %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the two epipoles
%---------------------%
[e1, e2] = e1_e2_from_F(F);

% compute the two 3-by-4 projection matrices
% result 9.15 Hartley, page 256
%
% note that for projective reconstruction we can choose the 3-by-4
% projection matrix P2 to have the form [ skew(e2)*F12+e2*alpha'  beta*e2 ], where
% alpha is an arbitrary 3-column vector and beta is an arbitrary scalar.  Here, for the rectification
% to work, the 3-by-3 submatrix of P2 (ie. the 1st three columns of P2) must be
% non-singular.  So, we adjust the vector alpha to make this submatrix non-singular.
% Since beta is also arbitrary, we can set beta to 1.
%
%---------------------%
P1    = [eye(3) zeros(3,1)];
alpha = .01*ones(1,3); % alpha is to make P2 non-singular
beta  = 1; % beta is an arbitrary scalar
P2    = [skew_r(e2)*F+e2*alpha  beta*e2];