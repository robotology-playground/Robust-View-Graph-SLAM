function a = R2a(R)
%function a = R2a(R)
%
% INPUT:  R - rotation matrix
% OUTPUT: a - Euler angle (radians)
%
% Adapted from the 'xyz' form in Matlab's dcm2angle.
%
% Tim Bailey 2012.


[r1 r2 r3] = threeaxisrot( -R(3,2,:), R(3,3,:), R(3,1,:), ...
                           -R(2,1,:), R(1,1,:));
a = [r1(:) r2(:) r3(:)]';

function [r1 r2 r3] = threeaxisrot(r11, r12, r21, r31, r32)
% find angles for rotations about X, Y, and Z axes
r1 = atan2( r11, r12 );
r2 = asin( r21 );
r3 = atan2( r31, r32 );
