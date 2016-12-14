% UNDISTPOINTINV - undistorts the point by the inverse mapping
%
% Usage: 
%           [x1, y1] = UndistPointInv(x2, y2, k1, k2, k3, p1, p2)
%
%
% Input:
%       x2, y2 : (normalized) distorted points
%       k1 : radial distortion coefficient
%       k2 : radial distortion coefficient
%       k3 : radial distortion coefficient
%       p1 : tangential distortion coefficient
%       p2 : tangential distortion coefficient
%
% Output:
%       x1, y1 : (normalized) undistorted points
%
%
% cf. Lens Distortion Model:
%       x2 = x1.*(1+k1.*r.^2 + k2.*r.^4 + k3.*r.^6) + 2.*p1.*x1.*y1 + p2.*(r.^2 + 2.*x1.^2)
%       y2 = y1.*(1+k1.*r.^2 + k2.*r.^4 + k3.*r.^6) + p1.*(r.^2 + 2.*y1.^2) + 2.*p2.*x1.*y1
%        where r.^2 = x1.^2 + y1.^2
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% May 2011 - Original version.


function [x1, y1] = UndistPointInv(x2, y2, k1, k2, k3, p1, p2)

% Initial Point
x1 = x2;
y1 = y2;

% Newton Method
IterMax = 20;
ReErrTh = 0.0001;   % Relative Error Threshold
TolTh   = 0.0001;   % Tolerence Threshold
for n = 1:IterMax
    r = sqrt(x1.^2 + y1.^2);

    % Lens Distortion Model
    f1 = x1.*(1+k1.*r.^2 + k2.*r.^4 + k3.*r.^6) + 2.*p1.*x1.*y1 + p2.*(r.^2 + 2.*x1.^2) - x2;
    f2 = y1.*(1+k1.*r.^2 + k2.*r.^4 + k3.*r.^6) + p1.*(r.^2 + 2.*y1.^2) + 2.*p2.*x1.*y1 - y2;

    % Partial Derivatives
    f1x = k1.*(x1.^2 + y1.^2) + k2.*(x1.^2 + y1.^2).^2 + k3.*(x1.^2 + y1.^2).^3 + ...
        6.*p2.*x1 + 2.*p1.*y1 + x1.*(2.*k1.*x1 + 4.*k2.*x1.*(x1.^2 + y1.^2) + 6.*k3.*x1.*(x1.^2 + y1.^2).^2) + 1;
    f1y = 2.*p1.*x1 + 2.*p2.*y1 + x1.*(2.*k1.*y1 + 4.*k2.*y1.*(x1.^2 + y1.^2) + 6.*k3.*y1.*(x1.^2 + y1.^2).^2);

    f2x = 2.*p1.*x1 + 2.*p2.*y1 + y1.*(2.*k1.*x1 + 4.*k2.*x1.*(x1.^2 + y1.^2) + 6.*k3.*x1.*(x1.^2 + y1.^2).^2);
    f2y = k1.*(x1.^2 + y1.^2) + k2.*(x1.^2 + y1.^2).^2 + k3.*(x1.^2 + y1.^2).^3 + ...
        2.*p2.*x1 + 6.*p1.*y1 + y1.*(2.*k1.*y1 + 4.*k2.*y1.*(x1.^2 + y1.^2) + 6.*k3.*y1.*(x1.^2 + y1.^2).^2) + 1;

    % Jacobian Matrix
    J = [f1x f1y; f2x f2y];
    F = [f1; f2];

    delta = -(J\F);

    x1_delta = x1 + delta(1);
    y1_delta = y1 + delta(2);

    ReErrX = abs((x1_delta - x1)./x1);
    ReErrY = abs((y1_delta - y1)./y1);

    TolX = abs(f1);
    TolY = abs(f2);

    x1 = x1_delta;
    y1 = y1_delta;
        
    if ((ReErrX < ReErrTh) & (ReErrY < ReErrTh)) | ...
            ((TolX < TolTh) & (TolY < TolTh))
        break;
    end
end