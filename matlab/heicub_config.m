function [options,encoders,floatingbase]=heicub_config()
%[options,encoders,floatingbase]=heicub_config()
%
%	For iCub@heidelberg.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014

% Restore original <Default>Properties of root,
% load default PATH, run STARTUP.m:
matlabrc;

% Toolbox functions
addpath('/home/tariq/Dev/mexopencv'); % mexopencv
addpath('/home/tariq/Dev/vlfeat-0.9.20/toolbox'); vl_setup(); % vlfeat
addpath('./batch/common'); %rvgslam
addpath('./batch/estimator_relative'); %rvgslam
addpath('./common'); %rvgslam
addpath('./mex'); %rvgslam
addpath('./test'); %rvgslam

% Robot related functions
addpath('./heicub'); 

DATA_DIR='/home/tariq/Documents/data/heicub/data_set1';
SAVE_DIR=strcat(DATA_DIR,'/run_',datestr(now,'yyyymmdd'));
CALB_DIR='/home/tariq/Documents/data/heicub/calib_20160913/img/stereo';

options=set_params(); % loads basic defaults
options=set_params(options,'folder',		DATA_DIR); % where is the data?
options=set_params(options,'save',			SAVE_DIR); % where to save the results?
options=set_params(options,'calib',			CALB_DIR); % where is the calibration file?
options=set_params(options,'freq',			10	); % frequency of acquisition for synchronisation
options=set_params(options,'first_image',	50	); % first image frame
options=set_params(options,'last_image',	60	); % last image frame
options=set_params(options,'steps',			2	); % frames resampling frequency
options=set_params(options,'verbose',		0	); % show verbose during data acquisition
 
[options,encoders,floatingbase]=set_images(options); % data acquisition and synchronisation
 
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
 
options=set_params(options,'optimiser'); % loads the pwg_optimiser defaults
options=set_params(options,'ncams',			5	); % PWGOPTIMISER: number of cams in each bundle
options=set_params(options,'nkeys',			1	); % PWGOPTIMISER: number of keyframes in each bundle
options=set_params(options,'nview',			10	); % PWGOPTIMISER: minimum number of views to accept a 3D point
options=set_params(options,'sigma_r',		1.0	); % PWGOPTIMISER: image measurements noise
options=set_params(options,'gateratio',		0.2	); % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
options=set_params(options,'verbose',		1	); % PWGOPTIMISER: show verbose during optimisation
options=set_params(options,'maxitr',		100	); % PWGOPTIMISER: maximum number of iterations in case no convergence
