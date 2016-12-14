function [H, e] = test_jacobian(xf, xc)

addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/matlab/utilities');

% select 1 point to evaluate
numpts = size(xf, 2); n = randperm(numpts);

% camera parameters
[K_left,K_right,d_left,d_right,Rc,tc]=test_param( );

% analytical jacobian
H = observation_model_jacobian(xc, xf(:,n(1)), K_left);

% numerical jacobian
J1 = numerical_jacobian_i(@observation_model, [], 1, [], xc, xf(:,n(1)), K_left);
J2 = numerical_jacobian_i(@observation_model, [], 2, [], xc, xf(:,n(1)), K_left);
J = [J1 J2];

% error metric
e = max(max(abs((H-J))));

% plots
%figure;
%spy(abs(H-J)>1e-06);