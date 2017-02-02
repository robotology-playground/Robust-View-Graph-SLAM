function [options,encoders,floatingbase]=r1_config()
%[options,encoders,floatingbase]=r1_config()
%
%	For r1@iit.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014, updated 2016

% Restore original <Default>Properties of root,
% load default PATH, run STARTUP.m:
matlabrc;

% Toolbox functions
addpath('/home/tariq/Dev/mexopencv'); % mexopencv for features extraction
addpath('/home/tariq/Dev/vlfeat-0.9.20/toolbox'); vl_setup(); % vlfeat for matching
addpath('../../suitesparse/CHOLMOD/MATLAB'); % CHOLMOD, for spinv
addpath('./batch/common'); % rvgslam, estimation common functions
addpath('./batch/estimator_relative'); % rvgslam for relative pose estimation
addpath('./batch/estimator_global'); % rvgslam for global pose estimation
addpath('./common'); % rvgslam, linear algebra, robot localisation, and vision related functions
addpath('./mex'); % rvgslam, mex wrappers to ../optimise and ../graph C++ classes
addpath('./test'); % rvgslam, a collection of test functions, some are needed here.

% Robot related functions
addpath('./r1'); % robot related functions

% Data related paths
DATA_DIR='/path/to/data';
SAVE_DIR=strcat(DATA_DIR,'/run_',datestr(now,'yyyymmdd')); % a folder with date suffix is created
CALB_DIR='./icub/calib_xxxxxxxx'; % path to intrinsic and extrinsic calibration files

% Simulation parameters:
options=set_params(); % loads basic defaults
% data collection parameters
options=set_params(options,'folder',		DATA_DIR); % where is the data?
options=set_params(options,'save',			SAVE_DIR); % where to save the results?
options=set_params(options,'calib',			CALB_DIR); % where is the calibration file?
options=set_params(options,'freq',			10	); 	% frequency of acquisition for synchronisation
													% Very large frequency means less accurate left to right synchronisation
													% Very small frequency means more frame drops
options=set_params(options,'first_image',	101	); % where to start reading the acquisition
options=set_params(options,'last_image',	200	); 	% where to stop reading the acquisition
													% both first_image and last_image represent (stereo) pair numbers
options=set_params(options,'steps',			2	); % frames resampling frequency (next_frame = current_frame + steps)
options=set_params(options,'verbose',		0	); 	% show verbose during data acquisition
													% 0 - no verbose, 1 - text logging, 2 - plot data
[options,encoders,floatingbase]=set_images(options); % run data acquisition and synchronisation
													% read inside notes to understand the purpose of set_images
% vision parameters
options=set_params(options,'vision'); % loads the vision defaults
options=set_params(options,'kazethreshold',	0.001); % KAZE: edge detection threshold
options=set_params(options,'kazeratio',		10);	% KAZE: matching distance between best two matches?
options=set_params(options,'mindisp',		2.0	); % minimum disparity to accept a correspondence
options=set_params(options,'ransac',		200	); % RANSAC: number of iterations
options=set_params(options,'RANSAC_pixtol',	0.5	); % RANSAC: pixel tolerance
options=set_params(options,'mincorrnr',		50	); % minimum number of corners to include and edge
options=set_params(options,'mininlnr',		25	); % minimum number of outliers to compute a geometry
options=set_params(options,'gridmargin',	5.0	); % discarded image margins during features extraction
options=set_params(options,'gridhorizon',	5.0	); % discarded image horison during features extraction
% optimiser parameters
options=set_params(options,'optimiser'); % loads the pwg_optimiser defaults
options=set_params(options,'ncams',			20	); % PWGOPTIMISER: number of cams in each bundle
%options=set_params(options,'nkeys',			10	); % PWGOPTIMISER: number of keyframes in each bundle
options=set_params(options,'nview',			10	); % PWGOPTIMISER: minimum number of views to accept a 3D point
options=set_params(options,'sigma_r',		1.0	); % PWGOPTIMISER: image measurements noise
options=set_params(options,'gateratio',		0.3	); % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
options=set_params(options,'verbose',		1	); % PWGOPTIMISER: show verbose during optimisation
options=set_params(options,'maxitr',		100	); % PWGOPTIMISER: maximum number of iterations in case no convergence
