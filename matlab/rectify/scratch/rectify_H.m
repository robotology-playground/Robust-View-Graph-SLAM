function [H,newe,G,R] = rectification_H(e, orgin)

% computes the rectification transformation for the image using its epipole
%
%   [H,newe] = rectification_H(e)
%   [H,newe,G,R] = rectification_H(e)
%   [H,newe,G,R] = rectification_H(e,origin) computes the rectification
%   transformation, H, for the image which has the epipole given by
%   the 3-by-1 vector, e.  The output argument, H, maps the epipole
%   e to a point at infinity -- typically, [1; 0; 0] for corresponding
%   horizontal scanlines or [0; 1; 0] for corresponding vertical scanlines.
%
%   The two optional output arguments, G and R, satisfy the equation
%      H = G*R*T,
%   with
%   - R being the 3-by-3 2D rotation matrix that brings the epipole e to
%     align with the x- or y-axis;
%   - G being the transformation that brings the epipole e, after rotation
%     R is applied, to infinity.
%   - T is the transformation that would transfer the origin of the image
%     coordinate system to the last input argument, origin.
%
%   The implementation here follows that described in "Theory and Practice
%   of Projective Rectification" by R. I. Hartley, IJCV'99.
%

if nargin == 1
    % assume that the origin of the image coordinate system has been
    % set at the centre of the image buffer
    T = eye(3);
else
    % assume that the origin was at the top-left corner.  We now move
    % it to the centre of the image, whose size is given by imsize.
    T = [1 0 -origin(1); 0 1 -origin(2); 0 0 1];
end

[ang,newe,R] = rectify_angle(e);
Re = R*e;
if sum(abs(newe-[1;0;0])) == 0
    % bring epipole e to parallel to [1;0;0]
    G = [1 0 0; 0 1 0; -Re(3)/Re(1) 0 1];
else
    % bring epipole e to parallel to [0;1;0]
    G = [1 0 0; 0 1 0; 0 -Re(3)/Re(2) 1];
end
H = G*R*T;

return


