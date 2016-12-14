% function test
clear all; clc; close all;
cd('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/icub');
run('./set_folders');
save_folder='/home/tabuhashim/Documents/MATLAB/koroibot/stereo/data/blue/run_20150927';
load(strcat(save_folder,'/options'));

% settings
load_kin=1; % load pwg from kinematics
load_vis=1; % load pwg from estimator
recompute_graph=0; % reduce the cost by removing some constraints
new_thickness=20; % can not go higher than opt.thickness
                  % as no matching data will be available

% recompute relatives from kinematics? if needed
if load_kin
    load(strcat(save_folder,'/matches'));
    load(strcat(save_folder,'/fwd_kinematics'));
else
    load(strcat(save_folder,'/encoders'));
    load(strcat(save_folder,'/matches'));
    encoders.eyes=eyes;
    encoders.neck=neck;
    encoders.waist=waist;
    if recompute_graph;% modify camera graph
        param.thickness=min(new_thickness,opt.thickness);
        n=2*size(opt.cam_left.image,2);
        mask=batch_initialise_camera_graph(param,n);
        G=mask.*G;
    end
    %[pwg,P]=batch_constraints_from_kinematics(encoders,pwm,G);
    [C,P]=batch_constraints_from_kinematics_new(encoders,C,options);
end

% accomulation of relatives from kinematics
%N=2*size(options.cam_left.image,2);
%[pwg,G]=C_to_pwg(C,N);
%[qk,tk]=pose_generate_sequential(G,pwg,{'Red','Blue','Cyan'},1);
%[qk,tk]=pose_generate_sequential_new(C,{'Red','Green','Lime'},0);

% recompute relatives from vision? if needed
if load_vis
    load_vis;load(strcat(save_folder,'/pairwise_geom'));
else
    load(strcat(save_folder,'/kpts'));
    load(strcat(save_folder,'/fwd_kinematics'));
    % modify camera graph
    if recompute_graph;
        param.thickness=min(new_thickness,opt.thickness);
        n=2*size(opt.cam_left.image,2);
        mask=batch_initialise_camera_graph(param,n);
        G=mask.*G;
    end
    % vision options
    opt.mindisp=2;
    opt.minbase=10; % in mm
    opt.RANSAC_pixtol=2;
    opt.ransac=50;
    opt.bucketsize=[150,150];
    % estimator options
    opt.iteration=2;
    opt.gateinnov=chi_square_bound(.99,2);
    opt.gateresid=chi_square_bound(.90,2);
    opt.gatetrust=chi_square_bound(.90,2);
    opt.gateratio=0.3;
    opt.checkrank=true;
    opt.verbose=1; % verbose level: 0,1,2
    % run the estimator
    [pwg,G]=batch_Emat_v2(pwg,G,kpts,opt);
end



[y,Y,sw,x,C,Ct]=pose_graph_estimation_new(C);




% accomulation of relatives from vision
%N=2*size(options.cam_left.image,2);
%[pwg,G]=C_to_pwg(C,N);
%[qv,tv]=pose_generate_sequential(G,pwg,{'Red','Green','Lime'},1);
%[qv,tv]=pose_generate_sequential_new(C,{'Red','Green','Lime'},1);

% spanning tree
%[q,t,sptree]=demo_spanning_tree(G,pwg,opt);

% motion averaging
%N=2*size(options.cam_left.image,2);
%pwg=C_to_pwg(C,N);
%[q,t,sptree]=batch_motion(G,pwg,options);

%C=pwg_to_C(pwg);
%[q,t]=batch_motion_new(C,options);

% least-squares
%[q,t]=pose_graph_estimation_new(qs,ts,C);

%x=[q t];
%figure;
%show_poses(x,{'Red','Green'});
%show_entire_graph_new(C,x,'Lime');