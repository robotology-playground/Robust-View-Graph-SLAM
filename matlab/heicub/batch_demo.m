%
%/*
% * Copyright (C) 2016 iCub Facility - Istituto Italiano di Tecnologia
% * Author: Tariq Abuhashim
% * Date: Nov 2016
% * email: t.abuhashim@gmail.com
% * Acknowledgement: This research has received funding from the European Union’s
% * Seventh Framework Programme for research, technological development and demonstration under
% * grant agreement No. 611909(KoroiBot).
% * CopyPolicy: Released under the terms of the LGPLv2.1 or later, see LGPL.TXT
% */
% 

addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/batch/common');
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/common');
addpath('/home/tariq/Documents/Robust-View-Graph-SLAM/matlab/test');
addpath('/home/tariq/Dev/mexopencv');
addpath('/home/tariq/Dev/vlfeat-0.9.20/toolbox'); vl_setup();

clc; clear all;

% 
[options, encoders, floatingbase] = heicub_config() ;

% Compute forward kinematics
Pkin = cameras_from_kinematics(encoders, floatingbase);

% Perform matching (or tracking) of image correspondences
[C,kpts] = build_camera_graph(options);

% Utilise kinematics or epipolar geometry to initialise constraints
C = initialise_graph_constraints(C,kpts,Pkin,options); 

% Refined sets of pair-wise geometry constraints using image correspondences
%C = optimise_pwg_constraints(C) ;

% Assign constraint weights


% Optimise for global camera poses using all the constraints in the graph
%C = optimise_graph_constraints(C) ;


