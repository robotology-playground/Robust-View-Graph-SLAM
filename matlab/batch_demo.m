% function batch_demo(robot)
%function batch_demo(robot)
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
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016
%
% if you receive and error: libstdc++.so.6: version `GLIBCXX_x.x.xx' not found:
% exit MATLAB
% rm -r /usr/local/MATLAB/R2016b/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6
% sudo ln -s /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/MATLAB/R2016b/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6
% restart MATLAB

% Restore original <Default>Properties of root,
% load default PATH, run STARTUP.m:
%matlabrc;
restoredefaultpath ;
clc;

robot = 'icub';

% Configure a robot
[options,encoders,floatingbase] = config_robot(robot);

% Compute forward kinematics
Pkin = cameras_from_kinematics(encoders,floatingbase);

% Perform matching (or tracking) of image correspondences
C = build_camera_graph(options);

% Utilise kinematics or epipolar geometry to initialise constraints
C = initialise_graph_constraints(C,Pkin,options); 

% Refined sets of pair-wise geometry constraints using image correspondences
%C = optimise_pwg_constraints(C,options);

% Assign constraint weights
%C = assign_constraint_weight(C);

% Optimise for global camera poses using all the constraints in the graph
%C = optimise_graph_constraints(C);

