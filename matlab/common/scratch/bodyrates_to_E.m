% transformation from body rotation rates to euler rotation rates

function E = bodyrates_to_E(a)

% the matrix E, converts from body angular rates to euler angular rates,
% given previous euler angle states (a).
% a: should be in randians, and in an order of [roll pitch yaw], 'XYZ'
% E = [a1 a2 a3]
%     [b1 b2 b3]
%     [c1 c2 c3]
%
% a = [roll pitch yaw] = [r p q]

a1 = 1;
a2 = sin(a(1))*tan(a(2));
a3 = cos(a(1))*tan(a(2));
b1 = 0;
b2 = cos(a(1));
b3 = -sin(a(1));
c1 = 0;
c2 = sin(a(1))*sec(a(2));
c3 = cos(a(1))*sec(a(2));

E = [a1 a2 a3; b1 b2 b3; c1 c2 c3];