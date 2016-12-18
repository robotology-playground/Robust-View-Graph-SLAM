function [P, A0, H0, A_left, A_right] = cameras_from_kinematics(encoders, floatingbase)
% [P, A0, H0, A_left, A_right] = cameras_from_kinematics(encoders, floatingbase)
%
% this function takes inputs as the pairwise_geometry of set of images and
% opt structure. The function will modefy and output a new
% pairwise_geometry where projection matrices are replaced by the forward
% kinematics of the robot. It also outputs the new camera graph, camera
% projection matrices relative to a reference, and the reference projection
% matrix.
%
% INPUTS:
%--------
% encoders		angles in the order of eyes, neck, waist
% floatingbase	ZMP computed root poses (note: could be unreliable)
%
% OUTPUTS:
%--------
% P			camera matrices from forward kinematics
% A0		reference camera matrix
% H0		reference frame transformation
% A_left	left camera matrix to current root frame
% A_right	right camera matrix to current root frame
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
%
% Tariq Abuhashim
% started: September 2014
%
% iCub - Koroibot

eyes = encoders([1 2 3],:);
neck = encoders([4 5 6],:);
waist = encoders([7 8 9],:);

if size(neck, 2)>1;
    disp('  ');
    disp('Camera poses from kinematics.');
end

H=[
0 1 0 0; 
0 0 1 0; 
1 0 0 0; 
0 0 0 1];

% compute camera matrices from forwark kinematics
P=cell(1,2*size(eyes,2));
A_left=cell(1,size(eyes,2));
A_right=cell(1,size(eyes,2));

c0=[zeros(3,1);R2w(H(1:3,1:3))'];
for i=1:size(eyes,2)
    % compute body-head-eye rotations and translation from kinematics
    [A,B]=forward_kinematics(waist(1:3,i),neck(1:3,i),eyes(1:3,i));
    % copy transformations into a new cell
    A_left{i}=scale_to_meters(A);
    A_right{i}=scale_to_meters(B);
    if nargin>1;% project to floating-base coordinates as reference
		ci=floatingbase(:,i);
		% Left camera
		cj=camera_matrix_to_pose(A_left{i});
		cj=transform_to_global_w(cj,c0);
		cji=transform_to_global_w(cj,ci);
		P{2*i-1}=[w2R(cji(4:6)) cji(1:3)];
		% right camera
		cj=camera_matrix_to_pose(A_right{i});
		cj=transform_to_global_w(cj,c0);
		cji=transform_to_global_w(cj,ci);
		P{2*i-0}=[w2R(cji(4:6)) cji(1:3)];
    else
		P{2*i-1}=A_left{i}(1:3,1:4);  % output A_left w.r.t the first joint
		P{2*i}=A_right{i}(1:3,1:4);  % output A_right w.r.t the first joint
    end
end

% normalise wrt the first camera?
A0=P{1};
ci=camera_matrix_to_pose(A0); % using kinematics
for i=1:size(P,2)
    cj=camera_matrix_to_pose(P{i});
    cij=transform_to_relative_w(cj,ci);
    %P{i}=[w2R(cij(4:6)) cij(1:3)];
end

% plot ?
if 0
   plot_cameras_floating(P, floatingbase);
   %pause
end
%
%
function P=scale_to_meters(P)
P(1:3,4)=P(1:3,4)/1000; % scale to meters?

function x=camera_matrix_to_pose(P)
x=[P(1:3,4);R2w(P(1:3,1:3))'];

function plot_cameras_floating(Pkin,floatingbase)
% COM and eyes
% floating-base
figure; hold on;
plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:));
plot3(floatingbase(1,1),floatingbase(2,1),floatingbase(3,1),'bO');
axis equal;grid on;box on;
axis tight;
% camera from kinematics
A0=Pkin{1}; pose=zeros(6,size(Pkin,2));
%ci=camera_matrix_to_pose(Pkin{1}); % using kinematics
for i=1:size(Pkin,2)
    %cj=camera_matrix_to_pose(Pkin{i});
    %pose(:,i)=transform_to_relative_w(cj,ci); % convert to relative (z is forward)
    pose(:,i)=camera_matrix_to_pose(Pkin{i});
end
% camera position
%figure; hold on;
plot3(pose(1,1:2:end),pose(2,1:2:end),pose(3,1:2:end),'r+');
plot3(pose(1,2:2:end),pose(2,2:2:end),pose(3,2:2:end),'g+');
plot3(pose(1,1),pose(2,1),pose(3,1),'bO');
axis equal; grid on; box on;
axis tight;
% camera and floating-base orientation
figure; hold on;
subplot(3,1,1); hold on; plot(pose(4,1:2:end)*180/pi,'r'); plot(floatingbase(4,:)*180/pi);
subplot(3,1,2); hold on; plot(pose(5,1:2:end)*180/pi,'r'); plot(floatingbase(5,:)*180/pi);
subplot(3,1,3); hold on; plot(pose(6,1:2:end)*180/pi,'r'); plot(floatingbase(6,:)*180/pi); 
% stereo baseline (relatives)
baseline=zeros(6,size(pose,2)/2);
k=0;
for i=2:2:size(pose,2)
    k=k+1;
    baseline(:,k)=transform_to_relative_w(pose(:,i),pose(:,i-1));
end
figure;
subplot(3,1,1); plot(baseline(1,:));
subplot(3,1,2); plot(baseline(2,:));
subplot(3,1,3); plot(baseline(3,:));
figure;
subplot(3,1,1); plot(baseline(4,:));
subplot(3,1,2); plot(baseline(5,:));
subplot(3,1,3); plot(baseline(6,:));
