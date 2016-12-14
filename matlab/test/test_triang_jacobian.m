% this is a script to check the triangulation jacobian only
clc; clear all; close all;

Rr=eye(3); xr=[.068;0;0];
u1 = rand; v1 = rand; u2 = rand; v2 = rand;

% numerical jacobian
J1 = numerical_jacobian_i(@test_triangulate,[],3,[],Rr,xr,u1,v1,u2,v2);
J2 = numerical_jacobian_i(@test_triangulate,[],4,[],Rr,xr,u1,v1,u2,v2);
J3 = numerical_jacobian_i(@test_triangulate,[],5,[],Rr,xr,u1,v1,u2,v2);
J4 = numerical_jacobian_i(@test_triangulate,[],6,[],Rr,xr,u1,v1,u2,v2);
J = [J1 J2 J3 J4]

% analytical jacobian
H = test_triangulate_jacobian(Rr, xr, u1, v1, u2, v2);
full(H)

% comparison
max(max(abs((H-J))))
%figure;
%spy(abs(H-J)>1e-06);