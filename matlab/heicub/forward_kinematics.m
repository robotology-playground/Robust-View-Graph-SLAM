function [A_left, A_right, A_inertial, Chain, H_0] = forward_kinematics(waist, neck, eyes)

% Calculates the camera and inerital extrinsic matrices from robot's
% kinematics model.
%
% Inputs:
% -----
% waist : torso orientation in radians as [roll, pitch, yaw].
%           Example:
%                      waist = [0, 0, 0]*pi/180;
%                      waist(waist==0) = .00001;
% neck : head group orientation in radians as [roll, pitch, yaw].
%           Example:
%                      neck = [0, 0, 0]*pi/180;
%                      neck(neck==0) = .00001;
% eyes: eyes orientation in radians as [tilt, pan, vergence].
%           Example:
%                      eyes = [0, 0, 0]*pi/180;
%                      eyes(eyes==0) = .00001;
%
% Outputs:
% -------
% A_left : is the coordinate tranformation from waist to left eye (A = [R  t])
% A_right : is the coordinate tranformation from waist to right eye
% A_inertial : is the coordinate tranformation from waist to inertial
% Chain : contains individual coordinate transformations of the whole chain
% H_0 : body reference to first link, which is the waist pitch
%
% Notes:
% -----
% Link 0: torso pitch
% Link 1: torso roll
% Link 2: torso yaw
% Link 3: neck pitch
% Link 4: neck roll
% Link 5: neck yaw
% The convention of this file, follows this page:
%   http://eris.liralab.it/wiki/ICubHeadV2Kinematics
%   torso (pitch -> roll -> heading) -> neck (pitch -> roll -> heading) -> then eyes
% Which is different than the one on this page:
%   http://eris.liralab.it/wiki/ICubForwardKinematics
%   where the figure shows,   
%   torso (heading -> roll -> pitch) -> neck (pitch -> roll -> heading) -> then eyes
%
% Examples:
% --------
% To calculate the camera orientation given change in waist roll angle
%           DCM = R_12*R_23*R_34*R_45*R_56*R_67*R_78
% To calculate the camera orientation given change in neck roll angle
%           DCM = R_45*R_56*R_67*R_78
%
% Tariq Abuhashim
% started: October 2014
% t.abuhashim@gmail.com
%
% iCub - Koroibot

% torso
if nargin<1; waist=zeros(1, 3); end
wthet0 = waist(2); % torso pitch
wthet1 = waist(1); % torso roll
wthet2 = waist(3); % torso yaw
% neck
if nargin<2; neck=zeros(1, 3); end
theta0 = neck(2); % head pitch
theta1 = neck(1); % head roll
theta2 = neck(3); % head yaw
% eyes
if nargin<3; eyes=zeros(1, 3); end
theta3 = eyes(1); % tilt
theta4 = eyes(2) - eyes(3)/2; % (R = Vs - Vg/2)
theta5 = eyes(2) + eyes(3)/2; % (L = Vs + Vg/2)

% transformations of links 0 to 5
H_01 = dh_matrix(32,    0,      pi/2, wthet0 + 0) ;    % waist pitch to waist roll
H_12 = dh_matrix(0,    -5.5,    pi/2, wthet1 - pi/2);  % waist roll to waist yaw
H_23 = dh_matrix(0,    -223.3, -pi/2, wthet2 - pi/2);  % waist yaw to neck pitch
H_34 = dh_matrix(0,   0,      pi/2, theta0 + pi/2) ; % neck pitch to neck roll
H_45 = dh_matrix(0,     0,     -pi/2, theta1 - pi/2);  % neck roll to neck yaw
H_56 = dh_matrix(0, 0, -pi/2, theta2 + pi/2) ; % neck yaw to eye tilt

% right eye
H2_67 = dh_matrix(0,    34,    -pi/2, theta3 + 0);
H2_78 = dh_matrix(0,    0,      pi/2, theta4 - pi/2); % (x - right, y - down, z - out of image)
%H_89 = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];  % (x - out of image, y - right, z - down)
A_right = H_01*H_12*H_23*H_34*H_45*H_56*H2_67*H2_78;%*H_89;
%A_right = H2_78*H2_67*H_56*H_45*H_34*H_23*H_12*H_01;

% left eye 
H1_67 = dh_matrix(0,   -34,    -pi/2, theta3 + 0);
H1_78 = dh_matrix(0,    0,      pi/2, theta5 - pi/2); % (x - right, y - down, z - out of image)
%H_89 = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];  % (x - out of image, y - right, z - down)
A_left  = H_01*H_12*H_23*H_34*H_45*H_56*H1_67*H1_78;%*H_89;
%A_left  = H1_78*H1_67*H_56*H_45*H_34*H_23*H_12*H_01;

% inertial
H_56 = dh_matrix(18.5,  110.8, -pi/2, theta2 + pi/2);
H_6I = dh_matrix(0, 6.6, pi/2, 0);
A_inertial = H_01*H_12*H_23*H_34*H_45*H_56*H_6I;

% full chain of the left eye (stereo vision system reference frame)
Chain = cell(10, 1);
Chain{1} = H_01; % waist pitch to waist roll
Chain{2} = H_12; % waist roll to waist yaw
Chain{3} = H_23; % waist yaw to neck pitch
Chain{4} = H_34; % neck pitch to neck roll
Chain{5} = H_45; % neck roll to neck yaw
Chain{6} = H_56; % neck yaw to left and right eyes tilt 
Chain{7} = H1_67; % left eye tilt to left eye heading
Chain{8} = H1_78; % left eye heading to left eye reference
Chain{9} = H2_67; % left eye tilt to right eye heading
Chain{10} = H2_78; % left eye heading to right eye reference

% reference transformation to first link
H_0 = [0 -1 0 0; 0 0 -1 0; 1 0 0 0; 0 0 0 1]; % waist reference to waist pitch
