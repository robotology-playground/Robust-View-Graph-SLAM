%
%	For iCub@heidelberg.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014, updated 2016

% Robots :
fprintf(' - Using iCub@heidelberg.\n') ;
robot = 'heicub';

% '/path/to/data'
DATA_DIR = ' ' ; 

% '/path/to/robot/calib_date'
CALB_DIR = './heicub/calib_20160913';

% External functions :
MEXOPENCV_DIR = '/home/tariq/Dev/mexopencv' ; % '/path/to/mexopencv', for features extraction
VLFEAT_TOOLBOX_DIR = '/home/tariq/Dev/vlfeat/toolbox' ; % '/path/to/vlfeat/toolbox', for matching
MATLAB_CHOLMOD_DIR = '/home/tariq/Dev/SuiteSparse/CHOLMOD/MATLAB' ; % '/path/to/SuiteSparse/CHOLMOD/MATLAB', for spinv

% Simulation parameters :
FREQUENCY 		= 10 ; % frequency of acquisition for synchronisation
FIRST_IMAGE 	= 61 ; % where to start reading the acquisition
LAST_IMAGE 		= 100 ; % where to stop reading the acquisition
DOWN_SAMPLE 	= 2 ; % frames resampling frequency (next_frame = current_frame + DOWN_SAMPLE)
VIS_VERBOSE 	= 0 ; % show verbose during data acquisition ( 0-quiet, 1-text logging, 2-plot )

% Vision parameters
KAZE_THRESHOLD 			= 0.001 ; % KAZE: edge detection threshold
KAZE_RATIO 				= 10 ;	% KAZE: matching distance between best two matches?
MIN_DISPARITY 			= 2.0 ; % minimum disparity to accept a correspondence
RANSAC_ITR 				= 200 ; % RANSAC: number of iterations
RANSAC_PIX_TOL 			= 0.5 ; % RANSAC: pixel tolerance
MIN_MUM_OF_EDGE_POINTS 	= 50 ; % minimum number of corners to include and edge
MIN_MUM_OF_EDGE_INLIERS = 25 ; % minimum number of outliers to compute a geometry

% optimiser parameters
NUM_OF_CAMS 	= 4 ; % PWGOPTIMISER: number of cams in each bundle
NUM_OF_VIEWS 	= 10 ; % PWGOPTIMISER: minimum number of views to accept a 3D point
IMAGE_NOISE 	= 0.5 ; % PWGOPTIMISER: image measurements noise
GATE_RATIO 		= 0.2 ; % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
MAX_NUM_ITR 	= 100 ; % PWGOPTIMISER: maximum number of iterations in case no convergence
OPT_VERBOSE 	= 2 ; % PWGOPTIMISER: show verbose during optimisation ( 0-quiet, 1-text logging, 2-plot )
