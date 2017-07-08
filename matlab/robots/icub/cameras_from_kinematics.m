function [P, A0, H0, A_left, A_right] = cameras_from_kinematics(encoders, floatingbase)

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

eyes = encoders([1 2 3],:);
neck = encoders([4 5 6],:);
waist = encoders([7 8 9],:);

if size(neck, 2)>1;
    %disp('  ');
    disp(' - Computing camera poses from kinematics.');
end

% compute camera matrices from forwark kinematics
P=cell(1,2*size(eyes,2));
A_left=cell(1,size(eyes,2));
A_right=cell(1,size(eyes,2));
%S=eye(4); %S(4,4) = .001;  % can introduce scale
%H=[0 0 1 0;-1 0 0 0;0 1 0 0;0 0 0 1];
%H=[0 -1 0 0;0 0 1 0;1 0 0 0;0 0 0 1];

for i=1:size(eyes,2);
    % compute body-head-eye rotations and translation from kinematics
    [A,B,~,~,H0]=forward_kinematics(waist(1:3,i),neck(1:3,i),eyes(1:3,i));
    % copy transformations into a new cell
    A_left{i}=scale_to_meters(A);
    A_right{i}=scale_to_meters(B);
    if nargin>3;% project to floating-base coordinates as reference
        T=[w2R(floatingbase(4:6,i)),floatingbase(1:3,i);0 0 0 1];
        A=(T*H0)*A_left{i};%*S;
        P{2*i-1}=A(1:3,1:4);  % output A_left w.r.t the floating-base
        B=(T*H0)*A_right{i};%*S;
        P{2*i}=B(1:3,1:4);  % output A_right w.r.t the floating-base
    else
        P{2*i-1}=A_left{i}(1:3,1:4);  % output A_left w.r.t the first joint
        P{2*i}=A_right{i}(1:3,1:4);  % output A_right w.r.t the first joint
    end
end

% normalise?
A0=P{1};
% for i=1:size(P,2)
%     %C=[A0;0 0 0 1]\[P{i};0 0 0 1]; %warning('kinematics are in first camera coordinates frame');
%     C=[P{i};0 0 0 1]; %warning('kinematics are in waist coordinates frame');
%     P{i}=scale_to_meters(C(1:3,1:4));
% end

function P=scale_to_meters(P)
P(1:3,4)=P(1:3,4)/1000; % scale to meters?
