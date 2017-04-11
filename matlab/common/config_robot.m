function [options,encoders,floatingbase]=config_robot(ROBOT_NAME)
%[options,encoders,floatingbase]=config_robot(ROBOT_NAME)
%
%	Configure robots
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014, updated 2016

% setup for icub, heicub, and r1
assert(any(strcmp(ROBOT_NAME,{'r1','icub','heicub'})),['Unknown robot: ', ROBOT_NAME]);
switch ROBOT_NAME
	case 'icub'%iCub@iit
		icub_config ; % FIXME: not tested yet
	case 'heicub'%iCub@heidelberg
		heicub_config ;
	case 'r1'%r1@iit
		r1_config ; % FIXME: not tested yet
end

% Toolbox functions
addpath(MEXOPENCV_DIR); % mexopencv for features extraction
addpath(VLFEAT_TOOLBOX_DIR); vl_setup(); % vlfeat for matching
addpath(MATLAB_CHOLMOD_DIR); % CHOLMOD, for spinv
addpath('./batch/common'); % rvgslam, estimation common functions
addpath('./batch/estimator_relative'); % rvgslam for relative pose estimation
addpath('./batch/estimator_global'); % rvgslam for global pose estimation
addpath('./common'); % rvgslam, linear algebra, robot localisation, and vision related functions
addpath('./mex'); % rvgslam, mex wrappers to ../optimise and ../graph C++ classes
addpath('./test'); % rvgslam, a collection of test functions, some are needed here.

% Robot related functions
addpath(['./',ROBOT_NAME]); % robot related functions

% Data related paths
SAVE_DIR=strcat(DATA_DIR,'/run_',datestr(now,'yyyymmdd')); % a folder with date suffix is created

% Simulation parameters:
options=set_params(); % loads basic defaults
% data collection parameters
options=set_params(options,'folder', DATA_DIR); % where is the data?
options=set_params(options,'save', SAVE_DIR); % where to save the results?
options=set_params(options,'calib', CALB_DIR); % where is the calibration file?
options=set_params(options,'freq', FREQUENCY); % frequency of acquisition for synchronisation
				% Very large frequency means less accurate left to right synchronisation
				% Very small frequency means more frame drops
options=set_params(options,'first_image', FIRST_IMAGE); % where to start reading the acquisition
options=set_params(options,'last_image', LAST_IMAGE); % where to stop reading the acquisition
				% both first_image and last_image represent (stereo) pair numbers
options=set_params(options,'steps', DOWN_SAMPLE); % frames resampling frequency (next_frame = current_frame + steps)
options=set_params(options,'verbose', VIS_VERBOSE); % show verbose during data acquisition
				% 0 - no verbose, 1 - text logging, 2 - plot data
[options,encoders,floatingbase]=set_images(options); % run data acquisition and synchronisation
				% read inside notes to understand the purpose of set_images
% vision parameters
options=set_params(options,'vision'); % loads the vision defaults
options=set_params(options,'detector', DETECTION_METHOD);
options=set_params(options,'detector_param', DETECTION_PARAMS);
options=set_params(options,'mindisp', MIN_DISPARITY); % minimum disparity to accept a correspondence
%options=set_params(options,'ransac', RANSAC_ITR); % RANSAC: number of iterations
%options=set_params(options,'RANSAC_pixtol', RANSAC_PIX_TOL); % RANSAC: pixel tolerance
options=set_params(options,'mincorrnr', MIN_MUM_OF_EDGE_POINTS); % minimum number of corners to include and edge
options=set_params(options,'mininlnr', MIN_MUM_OF_EDGE_INLIERS); % minimum number of outliers to compute a geometry
options=set_params(options,'gridmargin', 5.0); % discarded image margins during features extraction
options=set_params(options,'gridhorizon', 5.0); % discarded image horison during features extraction
% optimiser parameters
options=set_params(options,'optimiser'); % loads the pwg_optimiser defaults
options=set_params(options,'ncams', NUM_OF_CAMS); % PWGOPTIMISER: number of cams in each bundle
%options=set_params(options,'nkeys', 10); % PWGOPTIMISER: number of keyframes in each bundle
options=set_params(options,'nview', NUM_OF_VIEWS); % PWGOPTIMISER: minimum number of views to accept a 3D point
options=set_params(options,'sigma_r', IMAGE_NOISE); % PWGOPTIMISER: image measurements noise
options=set_params(options,'gateratio', GATE_RATIO); % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
options=set_params(options,'maxitr', MAX_NUM_ITR); % PWGOPTIMISER: maximum number of iterations in case no convergence
options=set_params(options,'verbose', OPT_VERBOSE); % PWGOPTIMISER: show verbose during optimisation

% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
% Query the last warning to acquire the identifier.  For example: 
% warnStruct = warning('query', 'last');
% msgid_integerCat = warnStruct.identifier
% msgid_integerCat =
%    MATLAB:concatenation:integerInteraction
warning('off','Images:initSize:adjustingMag');
