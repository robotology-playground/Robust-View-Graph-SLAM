function [zen, azm] = euclideanToSphere(poses)
% Convert (x,y,z) OpenGL coord system points on the unit sphere S2 to
% zenith/azimuth pairs. Angles are given in a compass frame:
% zenith=0 is always straight up (y-axis). Azimuth=0 points down -z, and
% Azimuth=pi/2 points down +x.
zen = acos(poses(:,2));
azm = atan2(poses(:,1), -poses(:,3));
return