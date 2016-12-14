function [pg,Rg] = transform2global(pr,Rr, pb,Rb)
%function [pg,Rg] = transform2global(pr,Rr, pb,Rb)
%
% Let {pb,Rb} be a frame defined with respect to the global frame, and
% {pr,Rr} a frame defined with respect to {pb,Rb}. This function transforms
% {pr,Rr} into {pg,Rg} in the global frame.
%
% Tim Bailey 2012.

pg = Rb*pr + pb;
Rg = Rb*Rr;
