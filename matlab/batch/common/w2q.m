function q = w2q(w)
%
% converts from euler angles to quaternions
% By Tariq Abuhashim
%
wtype = 0;
if size(w,1) == 3;
    w = w'; wtype = 1;
    if size(w,2) == 3; 
        warning('undetermined Euler vector is inverted to become row vector');
    end
end
w = w/2;
cw1 = cos(w(:,1)); sw1 = sin(w(:,1));
cw2 = cos(w(:,2)); sw2 = sin(w(:,2));
cw3 = cos(w(:,3)); sw3 = sin(w(:,3));

q(:,1) =  cw1.*cw2.*cw3 +  sw1.*sw2.*sw3;
q(:,2) = -cw1.*sw2.*sw3 +  cw1.*cw2.*sw3;
q(:,3) =  cw1.*cw2.*sw3 +  sw1.*cw2.*sw3;
q(:,4) =  cw1.*cw2.*sw3 -  sw1.*cw2.*sw3;

q = qnorm(q);
if wtype == 1; q = q'; end
return
