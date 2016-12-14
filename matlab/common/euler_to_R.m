% DCM

function R = euler_to_R(a)

% converts from euler angles (a) to direction cosine matrix (R)
% a: should be in randians, and in an order of [roll pitch yaw], 'XYZ'
% R = [a1 a2 a3]
%     [b1 b2 b3]
%     [c1 c2 c3]
%
% a = [roll pitch yaw] = [r p q]

c = cos(a);
s = sin(a);

a1 = c(3)*c(2);
a2 = c(3)*s(2)*s(1) - s(3)*c(1);
a3 = c(3)*s(2)*c(1) + s(3)*s(1);
b1 = s(3)*c(2);
b2 = s(3)*s(2)*s(1) + c(3)*c(1);
b3 = s(3)*s(2)*c(1) - c(3)*s(1);
c1 =-s(2);
c2 = c(2)*s(1);
c3 = c(2)*c(1);

R = [a1 a2 a3; b1 b2 b3; c1 c2 c3]; % standard DCM