function [ang,newe,R] = rectification_angle(e)

% computers the rectificaton angle, ang, required to bring the epipole e
% to infinity.
%
%   [ang,newe,R] = rectificationAngle(e) also returns the following
%   output arguments:
%   - the new epipole, newe, which would either be [1;0;0] or [0;1;0]
%     depending on the position of the epipole e.
%   - the image plane (2D) rotation matrix, R, required for rectification_H.m
%

e = e / norm(e,2);
% if e(3) is very small then the epipole must be at infinity already
if abs(e(3)) < 1E-20
    ang = 0;
    if abs(e(1)) > abs(e(2))
        newe = [1;0;0];
    else
        newe = [0;1;0];
    end
else
    theta = atan2(e(2)/e(3), e(1)/e(3))*ones(5,1);
    % below are 5 different angles to be considered.  If theta is
    % close to any of the first three angles then the rectification
    % should align the horizontal scanlines (new epipole should be
    % [1;0;0].  If theta is close to the last two angles then
    % the epipolar lines should be vertical (new epipole should be
    % [0;1;0].
    angles = [0; -pi; pi; -pi/2; pi/2];
    diff = angles - theta;
    [minval,minidx] = min(abs(diff));
    if minidx <= 3
        newe = [1;0;0];
    else
        newe = [0;1;0];
    end
    
    % compute the rotation angle require to bring the epipole e
    % to the x- or y-axis of the image
    ang = diff(minidx);
end

if nargout > 2
    R = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
end

return
