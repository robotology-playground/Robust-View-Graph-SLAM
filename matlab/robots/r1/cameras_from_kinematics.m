function [P, A0, H0, A__left, A_right] = cameras_from_kinematics(options, floatingbase)

%--------
% this function takes inputs as the pairwise_geometry of set of images and
% opt structure. The function will modefy and output a new
% pairwise_geometry where projection matrices are replaced by the forward
% kinematics of the robot. It also outputs the new camera graph, camera
% projection matrices relative to a reference, and the reference projection
% matrix.
%
% INPUTS:
%--------
% eyes, neck, waist: angles from encoders
%
% OUTPUTS:
%--------
% P   : projection matrices from forward kinematics
% A0  :
% H0  :
% A_left
% A_right
%
% To test the relative localisation:
%--------
%
% a = [rand(3,4); 0 0 0 1]; % is the rot-tra between a and reference in
%           reference coordinate frame
% b = [rand(3,4); 0 0 0 1]; % is the rot-tra between b and reference in
%           reference coordinate frame
% Then, the rot-tran between b and a in a reference coordinate frame can be
% done in two ways:
%
%   1.  b = a\b,
%       a = eye(4);
%
%   2.  b = [a(1:3,1:3)\b(1:3,1:3) a(1:3,1:3)\(b(1:3,4)-a(1:3,4)); 0 0 0 1]
%       a = eye(4);
%--------

if size(floatingbase, 2)>1;
    disp('  ');
    disp('Camera poses from kinematics.');
end

% compute camera matrices from forwark kinematics
P=cell(1,2*size(floatingbase,2));
A__left=cell(1,size(floatingbase,2));
A_right=cell(1,size(floatingbase,2));

% Note that, for R1, there is a rigit transformation between the wheel-base
% and the stereo center:
%
% R = I; t = [.013; 0.0; .961];
%
% This value was computed using Gazebo
%
% reference transformation to first link
%H0 = [0 0 1 0; 1 0 0 0; 0 1 0 0; 0 0 0 1]'; % waist reference to waist pitch
H0 = dh_matrix(.013, .961, -90*pi/180, -90*pi/180);
% split baseline between left and right
w = R2w(options.R);
R1 = w2R(w/2); t1 = options.R'*options.T/2; % left
R2 = w2R(-w/2); t2 = options.R'*options.T/2; t2(1)=-t2(1); % right

for i = 1 : size(floatingbase,2)
    % compute body-head-eye rotations and translation from kinematics
    %[A,B,~,~,H0]=forward_kinematics(waist(1:3,i),neck(1:3,i),eyes(1:3,i));
    % copy transformations into a new cell
    A__left{i}=scale_to_meters([R1 t1; 0 0 0 1]); % assumed sitting on cer-base reference (only for now)
    A_right{i}=scale_to_meters([R2 t2; 0 0 0 1]);
    if nargin>1;% project to floating-base coordinates as reference
        T=[w2R([0 0 pi_to_pi(floatingbase(6,i))]), floatingbase(1:3,i); 0 0 0 1]; % heading, x, y
        P{2*i-1}=T*H0*A__left{i};  % output A_left w.r.t the floating-base
        P{2*i+0}=T*H0*A_right{i};  % output A_right w.r.t the floating-base
    else
        P{2*i-1}=A__left{i};  % output A_left w.r.t the first joint
        P{2*i+0}=A_right{i};  % output A_right w.r.t the first joint
    end
end

% normalise?
A0=P{1};
%for i = 1 : size(P,2)
%    P{i}=A0\P{i}; %warning('kinematics are in first camera coordinates frame');
%end

function P=scale_to_meters(P)
P(1:3,4)=P(1:3,4)/1000; % scale to meters?