function H = map_to_infinity( x, cx, cy )

% Given a point and the desired origin (point of minimum projective
% distortion), compute a homograph H = G*R*T taking the point to the
% origin, rotating it to align with the X axis, then mapping it to
% infinity.

% compute a translation homography T taking the origin to the center
T = [ 1 0 -cx;
      0 1 -cy;
      0 0  1];
x = T * x;

% find a rotation homography R that rotates the translated epipole x to the 
% X axis
% Now rotate the translated x to align with the X axis.
cur_angle = atan2( x(2), x(1) );
R = [ cos( -cur_angle ), -sin( -cur_angle ), 0
      sin( -cur_angle ),  cos( -cur_angle ), 0
      0,                  0, 1 ];
x = R * x;

% calculate the transform mapping the rotated and translated epipole x to 
% infinity
if abs( x(3)/norm(x) ) < 1e-6
    % It's already at infinity
    G = eye(3);
else
    f = x(1)/x(3);
    G = [ 1   0  0
          0   1  0
        -1/f  0  1 ];
end;

H = G*R*T;