function w = q2w(q)
%
% converts from quaternions to euler angles
% By Tariq Abuhashim
%
qtype = isq(q); 
if qtype == 0; error('input must be a quaternion or quaternion vector');
elseif qtype == 1; q = q';
elseif qtype == 2; % do nothing
elseif qtype == 3; error('component q`s are either columns or rows (indeterminate)');
end
if isnormq(q) ~= 2;
    q = qnorm(q);
end
w(:,1) = atan2(2*(q(:,1)*q(:,2)+q(:,3)*q(:,4)),(1-2*q(:,2).^2-2*q(:,3).^2) );
w(:,2) = asin( 2*(q(:,1)*q(:,3)-q(:,4)*q(:,2)));
w(:,3) = atan2(2*(q(:,1)*q(:,4)+q(:,2)*q(:,3)),(1-2*q(:,3).^2-2*q(:,4).^2) );
if qtype == 1; w = w'; end
return
