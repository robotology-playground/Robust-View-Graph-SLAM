%
% The two datasets data_02 and data_03 contain 1 meter of walking data
% each. This code combines both data sets into one single experiment.
%
% Tariq Abuhashim, 2016
%

clear; clc; close all; %dbstop if error;
cd ('/home/tabuhashim/Dropbox/code/matlab/icub');
run ('./set_folders');

% myCluster=parcluster('local');
% myCluster.NumWorkers=3;
% %saveAsProfile(myCluster,'local2');
% %gcp;

% walking experiment (2 meters of walking combining data_02 and data_03
DATA_DIR_1 = '~/Documents/data/WalkingDatasets_1/data_02';
DATA_DIR_2 = '~/Documents/data/WalkingDatasets_1/data_03';
SAVE_DIR = ['~/Documents/data/WalkingDatasets_1/results','/run_',datestr(now,'yyyymmdd')];
CALB_DIR{1} = '~/Documents/data/WalkingDatasets_1/data_calibration/img/left/results/run_1';
CALB_DIR{2} = '~/Documents/data/WalkingDatasets_1/data_calibration/img/right/results/run_1';

% options to verify all data logs
options = set_params(); % loads basic defaults
options = set_params(options, 'save', SAVE_DIR);
options = set_params(options, 'calib', CALB_DIR);
options = set_params(options, 'verbose', 0);

% first meter
options = set_params(options, 'folder', DATA_DIR_1);
options = set_params(options, 'freq', 30);
options = set_params(options, 'first_image', 861);%861
options = set_params(options, 'last_image', 5860);%5860
options = set_params(options, 'steps', 2); % set to 1 to check data limits first
[options_1, encoders_1, floatingbase_1] = set_images (options);

% second meter
options = set_params(options, 'folder', DATA_DIR_2);
options = set_params(options, 'freq', 30);
options = set_params(options, 'first_image', 540);
options = set_params(options, 'last_image', 4700);
options = set_params(options, 'steps', 2); % set to 1 to check data limits first
[options_2, encoders_2, floatingbase_2] = set_images (options);

% combine both data sets
m = size(floatingbase_1, 2);
n = size(floatingbase_2, 2);
ref_pose = floatingbase_1(:, m);
ref_pose(3:6) = zeros(4,1);
floatingbase_2 = transform_to_global_w(floatingbase_2, ref_pose);
floatingbase = [floatingbase_1, floatingbase_2];
encoders = [encoders_1, encoders_2];
eyes = encoders(1:3,:);
neck = encoders(4:6,:);
waist = encoders(7:9,:);

% combine image names
k = 0;
for i = 1 : size(options_1.cam_left.image,2)
    k = k+1;
    options.cam_left.image{k} = strcat(options_1.img_folder{1}, options_1.cam_left.image{i});
end
for i = 1 : size(options_2.cam_left.image,2)
    k = k+1;
    options.cam_left.image{k} = strcat(options_2.img_folder{1}, options_2.cam_left.image{i});
end

k = 0;
for i = 1 : size(options_1.cam_right.image,2)
    k = k + 1;
    options.cam_right.image{k} = strcat(options_1.img_folder{2}, options_1.cam_right.image{i});
end
for i = 1 : size(options_2.cam_right.image,2)
    k = k + 1;
    options.cam_right.image{k} = strcat(options_2.img_folder{2}, options_2.cam_right.image{i});
end

% set up the image processing and PwgOptimiser
options.folder = []; % Empty becasuse image names included the folders when they were combined above
options.img_folder{1} = [];
options.img_folder{2} = [];
options = set_params(options, 'vision'); % loads the vision defaults
options = set_params(options, 'fastthreshold',10);
options = set_params(options, 'mindisp',2.5);
options = set_params(options, 'ransac',500);
options = set_params(options, 'mincorrnr',50);
options = set_params(options, 'gridmargin',5);
options = set_params(options, 'gridhorizon',5);
options = set_params(options, 'optimiser'); % loads the pwg_optimiser defaults
%options = set_params(options,'verbose',0);
options = set_params(options, 'ncams',10);
options = set_params(options, 'nkeys',10);
options = set_params(options, 'nview',5);
options = set_params(options, 'sigma_r',0.5);
options = set_params(options, 'gateratio',0.5);
if isfield(options, 'save');
    save(strcat(options.save, '/options'), 'options');
    save(strcat(options.save, '/encoders'), 'eyes', 'neck', 'waist');
    save(strcat(options.save, '/floatingbase'), 'floatingbase');
end

% get kinematics
[Pkin, A0] = cameras_from_kinematics(eyes, neck, waist, floatingbase);

% Check kinematics and floating-base solution
if 0
    % COM and eyes
    pose = zeros(6, size(Pkin,2));
    for i = 1 : size(Pkin,2)
        pose(1:3,i) = Pkin{i}(1:3,4);
        pose(4:6,i) = R2w(Pkin{i}(1:3,1:3))';
    end
    figure; subplot(1,2,1); hold on;
    plot3(floatingbase(2,10:end-10), floatingbase(1,10:end-10), floatingbase(3,10:end-10));
    axis equal; axis([-.4 .4 0 2.5]); grid on; box on;
    subplot(1,2,2); hold on;
    plot3(pose(2,1:2:end), pose(1,1:2:end), pose(3,1:2:end), 'r+');
    plot3(pose(2,2:2:end), pose(1,2:2:end), pose(3,2:2:end), 'g+'); 
    axis equal; axis([-.4 .4 0 2.5]); grid on; box on;
    xlabel('Eastern (m)'); ylabel('Forward (m)');
    % encoders
    
    % baseline
    baseline = zeros(6, size(pose, 2)/2);
    k = 0;
    for i = 2 : 2 : size(pose, 2)
        k = k + 1;
        baseline(:,k) = transform_to_relative_w(pose(:, i), pose(:, i-1));
    end
    figure;
    subplot(3,1,1); plot(baseline(1,:));
    subplot(3,1,2); plot(baseline(2,:));
    subplot(3,1,3); plot(baseline(3,:));
    return
end

if 0 % for a forward walking video.
    for i = 1 : 1 : size(options.cam_left.image,2)
        start = tic;
        im1 = imread(options.cam_left.image{i});
        im2 = imread(options.cam_right.image{i});
        imshow([im1 im2])
        %T = get(gca,'tightinset');
        %set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
        %set(gcf, 'menubar', 'none');
        pause((options.steps/30)-toc(start))
    end
    return
end

if 0 % run sparse bundle estimation
    [C, scan] = batch_motion_and_map_inverse_depth_v2(options, Pkin);
    %delete(gcp);
    return
end