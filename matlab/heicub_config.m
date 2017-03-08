function [options,encoders,floatingbase]=heicub_config(DATA_DIR)
%[options,encoders,floatingbase]=heicub_config(DATA_DIR)
%
%	For iCub@heidelberg.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014, updated 2016

% Restore search path to its factory-installed state
restoredefaultpath;

% Dependencies
addpath('./mexopencv'); % mexopencv for features extraction
addpath('../vlfeat/toolbox'); vl_setup(); % vlfeat for matching
addpath('../SuiteSparse/CHOLMOD/MATLAB'); % CHOLMOD, for spinv

% RVGSLAM functions
addpath('./batch/common'); % rvgslam, estimation common functions
addpath('./batch/estimator_relative'); % rvgslam for relative pose estimation
addpath('./batch/estimator_global'); % rvgslam for global pose estimation
addpath('./common'); % rvgslam, linear algebra, robot localisation, and vision related functions
addpath('./mex'); % rvgslam, mex wrappers to ../optimise and ../graph C++ classes
addpath('./test'); % rvgslam, a collection of test functions, some are needed here.

% Robot related functions
addpath('./heicub'); % robot related functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data related paths
SAVE_DIR=strcat(DATA_DIR,'/run_',datestr(now,'yyyymmdd')); % a folder with date suffix is created
CALB_DIR=[pwd '/heicub/calib_20160913']; % path to intrinsic and extrinsic calibration files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters:
options=set_params(); % loads basic defaults
% data collection parameters
options=set_params(options,'folder',		DATA_DIR); % where is the data?
options=set_params(options,'save',			SAVE_DIR); % where to save the results?
options=set_params(options,'calib',			CALB_DIR); % where is the calibration file?
options=set_params(options,'freq',			10	); 	% frequency of acquisition for synchronisation
													% Very large frequency means less accurate left to right synchronisation
													% Very small frequency means more frame drops
options=set_params(options,'first_image',	61	); % where to start reading the acquisition
options=set_params(options,'last_image',	100	); 	% where to stop reading the acquisition
													% both first_image and last_image represent (stereo) pair numbers
options=set_params(options,'steps',			2	); % frames resampling frequency (next_frame = current_frame + steps)
options=set_params(options,'verbose',		0	); 	% show verbose during data acquisition
													% 0 - no verbose, 1 - text logging, 2 - plot data
													
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% Get the data										
[options,encoders,floatingbase]=set_images(options); % run data acquisition and synchronisation
													% read inside notes to understand the purpose of set_images
													
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vision parameters
options=set_params(options,'vision'); % loads the vision defaults
options=set_params(options,'kazethreshold',	0.001); % KAZE: edge detection threshold
options=set_params(options,'kazeratio',		10);	% KAZE: matching distance between best two matches?
options=set_params(options,'mindisp',		2.0	); % minimum disparity to accept a correspondence
options=set_params(options,'ransac',		200	); % RANSAC: number of iterations
options=set_params(options,'RANSAC_pixtol',	0.5	); % RANSAC: pixel tolerance
options=set_params(options,'mincorrnr',		50	); % minimum number of matches to compute two-view geometries
options=set_params(options,'mininlnr',		25	); % minimum number of inliers to trust two-view results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimiser parameters
options=set_params(options,'optimiser'); % loads the pwg_optimiser defaults
options=set_params(options,'ncams',			4	); % PWGOPTIMISER: number of cams in each bundle
%options=set_params(options,'nkeys',			10	); % PWGOPTIMISER: number of keyframes in each bundle
options=set_params(options,'nview',			10	); % PWGOPTIMISER: minimum number of views to accept a 3D point
options=set_params(options,'sigma_r',		0.5	); % PWGOPTIMISER: image measurements noise
options=set_params(options,'gateratio',		0.2	); % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
options=set_params(options,'verbose',		2	); % PWGOPTIMISER: show verbose during optimisation
options=set_params(options,'maxitr',		100	); % PWGOPTIMISER: maximum number of iterations in case no convergence
