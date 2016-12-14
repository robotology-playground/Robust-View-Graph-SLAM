function xr=relpose(x1,x2)

% translate
xr(1,1)=x2(1,1)-x1(1,1);
xr(2,1)=x2(2,1)-x1(2,1);
xr(3,1)=x2(3,1)-x1(3,1);
% rotate
R1=w2R(x1(4:6,1));
%R1=get_rotation_matrix(x1(4:6,1),'rodrigues');
xr(1:3,1)=R1'*xr(1:3,1);

% orientation using subtraction (use the analytical jacobian in
% relpose_jacobian_6DoF.m)
xr(4:6,1)=x2(4:6,1)-x1(4:6,1);
xr(5)=pi_to_pi(xr(5));

% orientation using rotation matrix (if using this, you have to use
% the numerical jacobian in relpose_jacobian_6DoF.m)
%R2=w2R(x2(4:6,1));
%Rr=R1'*R2;
%xr(4:6,1)=R2w(Rr);
%xr(5)=pi_to_pi(xr(5));