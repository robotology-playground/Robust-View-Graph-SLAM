clear;clc;close all;dbstop if error;
cd('/home/tabuhashim/Dropbox/code/matlab/icub');
run('./set_folders');

% myCluster=parcluster('local');
% myCluster.NumWorkers=3; 
% %saveAsProfile(myCluster,'local2');
% %gcp;

% walking experiment
DATA_DIR='~/Documents/data/WalkingDatasets/data_02';
SAVE_DIR='~/Documents/data/WalkingDatasets/results/run_2';
CALB_DIR{1}='~/Documents/data/WalkingDatasets/data_calibration/img/left/results/run_1';
CALB_DIR{2}='~/Documents/data/WalkingDatasets/data_calibration/img/right/results/run_1';
% FREQ=30;
% FIMAGE=1001;
% LIMAGE=1400;
% STEPS=5;
% options to verify all data logs
options=set_params(); %loads basic defaults
options=set_params(options,'folder',DATA_DIR);
options=set_params(options,'save',SAVE_DIR);
options=set_params(options,'calib',CALB_DIR);
options=set_params(options,'freq',30);
options=set_params(options,'verbose',0);
%set_images(options); pause; %checks all the data logs
% options to run the demo
options=set_params(options,'first_image',861);
options=set_params(options,'last_image',5860);
options=set_params(options,'steps',10);
[options,encoders,floatingbase]=set_images(options);

options=set_params(options,'vision');
options=set_params(options,'fastthreshold',10);
options=set_params(options,'mindisp',10);
options=set_params(options,'ransac',500);
options=set_params(options,'gridmargin',5);
options=set_params(options,'gridhorizon',5);
options=set_params(options,'optimiser');
%options=set_params(options,'verbose',0);
options=set_params(options,'ncams',30);
options=set_params(options,'nkeys',30);
options=set_params(options,'nview',10);
options=set_params(options,'sigma_r',0.5);
options=set_params(options,'gateratio',0.2);

% get kinematics
load(strcat(options.save,'/encoders'));
[Pkin,A0]=cameras_from_kinematics(eyes,neck,waist,floatingbase);
if 0% Check kinematics and floating-base solution
    for i=1:size(Pkin,2)
        pose(1:3,i)=Pkin{i}(1:3,4);
        pose(4:6,i)=R2w(Pkin{i}(1:3,1:3))';
    end
    clf;
    plot3(floatingbase(1,:),floatingbase(2,:),floatingbase(3,:));
    hold on;
    plot3(pose(1,1:2:end),pose(2,1:2:end),pose(3,1:2:end),'r+');
    plot3(pose(1,2:2:end),pose(2,2:2:end),pose(3,2:2:end),'g+'); 
    axis equal;
    figure; baseline=zeros(6,size(pose,2)/2);
    k=0;
    for i=2:2:size(pose,2)
        k=k+1;
        baseline(:,k)=transform_to_relative_w(pose(:,i),pose(:,i-1));
    end
    subplot(3,1,1);plot(baseline(1,:));
    subplot(3,1,2);plot(baseline(2,:));
    subplot(3,1,3);plot(baseline(3,:));
    return;
end

if 1% run sparse bundle estimation
    [C,scan]=batch_motion_and_map_inverse_depth_v2(options,Pkin);
    %delete(gcp);
    return;
end