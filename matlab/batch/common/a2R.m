function R = a2R(a)
%function R = a2R(a)
%
% INPUT:  a - Euler angle (radians)
% OUTPUT: R - rotation matrix
%
% Adapted from the 'xyz' form in angle2dcm, from the commercial Matlab toolbox.
%
% Tim Bailey 2012.

R = zeros(3,3,size(a,2));
ca = cos(a);
sa = sin(a);

%     [          cy*cz, sz*cx+sy*sx*cz, sz*sx-sy*cx*cz]
%     [         -cy*sz, cz*cx-sy*sx*sz, cz*sx+sy*cx*sz]
%     [             sy,         -cy*sx,          cy*cx]

R(1,1,:) = ca(2,:).*ca(3,:);
R(1,2,:) = sa(1,:).*sa(2,:).*ca(3,:) + ca(1,:).*sa(3,:);
R(1,3,:) = -ca(1,:).*sa(2,:).*ca(3,:) + sa(1,:).*sa(3,:);
R(2,1,:) = -ca(2,:).*sa(3,:);
R(2,2,:) = -sa(1,:).*sa(2,:).*sa(3,:) + ca(1,:).*ca(3,:);
R(2,3,:) = ca(1,:).*sa(2,:).*sa(3,:) + sa(1,:).*ca(3,:);        
R(3,1,:) = sa(2,:);
R(3,2,:) = -sa(1,:).*ca(2,:);
R(3,3,:) = ca(1,:).*ca(2,:);
