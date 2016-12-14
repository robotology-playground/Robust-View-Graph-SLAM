% UNDISTPOINTFOR - undistorts the point by the forward mapping
%
% Usage:
%           [x2, y2] = UndistPointFor(x1, y1, k1, k2, k3, p1, p2)
%
%
% Input:
%       x1, y1 : (normalized) source point
%       k1 : radial distortion coefficient
%       k2 : radial distortion coefficient
%       k3 : radial distortion coefficient
%       p1 : tangential distortion coefficient
%       p2 : tangential distortion coefficient
%
% Output:
%       x2, y2 : (normalized) destination point
%
%
% cf. Lens Distortion Model:
%       x2 = x1*(1+k1*r^2 + k2*r^4 + k3*r^6) + 2*p1*x1*y1 + p2*(r^2 + 2*x1^2)
%       y2 = y1*(1+k1*r^2 + k2*r^4 + k3*r^6) + p1*(r^2 + 2*y1^2) + 2*p2*x1*y1
%        where r^2 = x1^2 + y1^2
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Sep. 2010 - Original version.


function [x2, y2] = UndistPointFor(x1, y1, k1, k2, k3, p1, p2)

r = sqrt(x1.^2 + y1.^2);

% Undistorted Point
x2 = x1.*(1+k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1.*x1.*y1 + p2*(r.^2 + 2*x1.^2);
y2 = y1.*(1+k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*y1.^2) + 2*p2.*x1.*y1;