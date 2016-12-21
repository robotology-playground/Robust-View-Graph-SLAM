% function batch_demo()
%function batch_demo()
%
% Performs view-graph slam 
% Requirements:
%	1- Robot configuaration file: for examples, see
%		config_heicub(); % for iCub@heidelberg
%		config_icub(); % for iCub@iit
%		config_r1(); % for R1@iit
%	2- Vlfeat, which contains basic features matching tasks
%	3- mexopencv, which contains various Opencv wrappers.
%	4- rvgslam, which contains our robust view-graph slam implementation.
%
%	For all robots.
%
% Tariq Abuhashim - 2016
% t.abuhashim@gmail.com
%
% iCub - Koroibot

addpath('/home/tariq/Dev/mexopencv'); % mexopencv
addpath('/home/tariq/Dev/vlfeat-0.9.20/toolbox'); vl_setup(); % vlfeat
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/batch/common'); %rvgslam
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/batch/estimator_relative'); %rvgslam
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/batch/estimator_global'); %rvgslam
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/common'); %rvgslam
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/test'); %rvgslam
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/mex'); %rvgslam

clc; clear all; % FIXME: temporary thing, remove after using this as a function

% FIXME: setup for icub, heicub, and r1
[options, encoders, floatingbase] = heicub_config() ;

% Compute forward kinematics
Pkin = cameras_from_kinematics(encoders, floatingbase) ;

% Perform matching (or tracking) of image correspondences
[C,kpts] = build_camera_graph(options) ;

% Utilise kinematics or epipolar geometry to initialise constraints
C = initialise_graph_constraints(C,kpts,Pkin,options) ; 

% Refined sets of pair-wise geometry constraints using image correspondences
C = optimise_pwg_constraints(C,kpts,options) ;

% Assign constraint weights


% Optimise for global camera poses using all the constraints in the graph
%C = optimise_graph_constraints(C) ;


