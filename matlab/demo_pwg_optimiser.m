%function demo_pwg_optimiser
% A demo of robust pair-wise constraints optimiser
% Tariq Abuhashim, 2016.
% iCub - Koroibot

clear; clc; close all;
dbstop if error; % allows debug mode when error is declared

robot = 'heicub';
NEW_OPTIONS = 1;

if  NEW_OPTIONS
    
    %DATA_DIR='~/Documents/data/blue/20150515_1'];
    %DATA_DIR='~/Documents/data/blue/20150901_3'];
    %CALB_DIR='~/Documents/data/blue/calib_20150422';
    %SAVE_DIR = ['~/Documents/data/blue','/run_',datestr(now,'yyyymmdd')] ;
    
    %DATA_DIR='~/Documents/data/blue/20151207_desk3';
    %%DATA_DIR='~/Documents/data/blue/20151207_moving_2';
    %CALB_DIR{1}='~/Documents/data/blue/calib_20151207/left/results/run_2';
    %CALB_DIR{2}='~/Documents/data/blue/calib_20151207/right/results/run_2';
    %SAVE_DIR = ['~/Documents/data/blue','/run_',datestr(now,'yyyymmdd')];
    
    if strcmp(robot, 'icub')
        cd('/home/tabuhashim/Dropbox/code/matlab/icub');
        run('./set_folders'); % includes other functions
        
        % insert icub config here
        [options,encoders,floatingbase]=icub_config;
        
        % get globals from kinematics
        options.root=0; % robot root, P{1} = A0;
        load(strcat(options.save,'/encoders')); % loads eyes, neck, waist and floatingbase
        [Pkin,A0]=cameras_from_kinematics(eyes,neck,waist,floatingbase);
        
    elseif strcmp(robot, 'r1')
        cd('/home/tariq/Documents/MATLAB/r1');
        run('./set_folders'); % includes other functions
        
        % insert r1 config here
        [options,floatingbase]=r1_config();
        
        % get globals from kinematics
        floatingbase=transform_to_relative_w(floatingbase,floatingbase(:,1));
        [Pkin,A0]=cameras_from_kinematics(options,floatingbase);
        
    elseif strcmp(robot, 'heicub')
        cd('/home/tariq/Documents/MATLAB/heicub');
        run('./set_folders'); % includes other functions
        
        % insert heicub config here
        [options,waist,floatingbase]=heicub_config();
        
        % get globals from kinematics
        % for pose lists only
        %Rbt = w2R(floatingbase(4:6,1));
        %for i = 1:size(floatingbase, 2)
        %    Rp = w2R(floatingbase(4:6,i));
        %    Rp = Rbt*Rp;
        %    floatingbase(4:6,i) = R2w(Rp)';
        %end
        %floatingbase(5,:) = pi_to_pi(floatingbase(5,:));
        %floatingbase=transform_to_relative_w(floatingbase,floatingbase(:,1));
        [Pkin,A0]=cameras_from_kinematics(waist,floatingbase);
        
    else
        error('Invalid robot ...')
        
    end
    
else
    
    if strcmp(robot, 'icub')
        %folder='/home/tabuhashim/Documents/data/blue/run_20160213';
        %folder='/home/tabuhashim/Documents/data/blue/run_20160224';
        %folder='/home/tabuhashim/Documents/data/blue/run_20160301';
        folder='/home/tabuhashim/Documents/data/blue/run_20160319';
        %folder='/home/tabuhashim/Documents/data/blue/run_20160327';
        load(strcat(folder,'/options'));
        
    elseif strcmp(robot, 'r1')
        error('No tests on r1 yet.')
        
    elseif strcmp(robot, 'heicub')
        error('No tests on heicub yet.')
        
    else
        error('Invalid robot ...')
        
    end
    
end

% Check kinematics and floating-base solution
if 0
    % COM and eyes
    pose=zeros(6,size(Pkin,2));
    for i = 1 : size(Pkin,2)
        pose(1:3,i)=Pkin{i}(1:3,4);
        pose(4:6,i)=R2w(Pkin{i}(1:3,1:3))';
    end
    
    figure; hold on;
    plot3(floatingbase(1,10:end-10),floatingbase(2,10:end-10),floatingbase(3,10:end-10));
    plot3(floatingbase(1,1),floatingbase(2,1),floatingbase(3,1),'bO');
    axis equal;grid on;box on;
    axis tight;
    %axis([-.4 .4 0 2.5]);
    xlabel('Eastern (m)');ylabel('Forward (m)');
    
    figure; hold on;
    plot3(pose(1,1:2:end),pose(2,1:2:end),pose(3,1:2:end),'r+');
    plot3(pose(1,2:2:end),pose(2,2:2:end),pose(3,2:2:end),'g+');
    plot3(pose(1,1),pose(2,1),pose(3,1),'bO');
    %plot3(-pose(2,1:2:end),-pose(3,1:2:end),pose(1,1:2:end),'r+');
    %plot3(-pose(2,2:2:end),-pose(3,2:2:end),pose(1,2:2:end),'g+');
    axis equal; grid on; box on;
    axis tight;
    %axis([-.4 .4 0 2.5]);
    xlabel('Eastern (m)');ylabel('Forward (m)');
    
    % baseline
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
    return
end

if 1 % M-views implementation
    options.xest=zeros(6,1); % one of the termination criteria
    [C,scan]=batch_motion_and_map_inverse_depth_v2(options,Pkin);
end

% 2-views implementation
%[kpts, desc] = batch_features(options);
%[C, G] = batch_matching_new(kpts, desc, options);
%[C, P] = batch_constraints_from_kinematics_new(encoders, C, options);
%[C, scan] = batch_Emat_v2_new(C, kpts, options, P);
