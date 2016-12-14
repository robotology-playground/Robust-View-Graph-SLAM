function [pr,Rr] = transform2relative(p1,R1, p2,R2)
%function [pr,Rr] = transform2relative(p1,R1, p2,R2)
%
% Let {p1,R1} and {p2,R2} be frames defined with respect to the global
% frame. This function transforms {p1,R1} into {pr,Rr} in the frame of
% {p2,R2}. 
%
% Tim Bailey, 2012.

R2t = inv(R2); % rotation matrices are orthogonal; transpose equals inverse
pr = R2t*(p1 - p2);
Rr = R2t*R1;
