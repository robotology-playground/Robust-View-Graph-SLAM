% DCM

function a = R_to_euler(R)

% converts from direction ce matrix (R) to euler angles (a).
% a: should be in randians, and in an order of [roll pitch yaw], 'XYZ'
% R = [a1 a2 a3]
%     [b1 b2 b3]
%     [c1 c2 c3]
%
% a = [roll pitch yaw] = [r p q]

a(1,1) = atan2(R(3,2), R(3,3));
a(2,1) = asin(-R(3,1));
a(3,1) = atan2(R(2,1), R(1,1));