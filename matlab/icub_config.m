%
%	For iCub@iit.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2014, updated 2016

% Robots :
fprintf(' - Using iCub@iit.\n') ;
robot = 'icub';

% '/path/to/data'
DATA_DIR = '/home/tariq/Documents/data/blue/20151207_desk3';

% '/path/to/robot/calib_date'
CALB_DIR = './icub/blue/calib_20151207';

% External functions :
% in order to compile mexopencv without installing opencv:
% expert ${OPENCV_DIR}/lib/pkgconfig into PKG_CONFIG_PATH in your .bashrc 
% export PKG_CONFIG_PATH=${OPENCV_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH
% to make all, use: 	make all MATLABDIR=/usr/local/MATLAB/R2017a
% to make contrib, use: make contrib MATLABDIR=/usr/local/MATLAB/R2017a
% to make clean, use: 	make clean MATLABDIR=/usr/local/MATLAB/R2017a
% If you get an error: Invalid MEX-file '/mexopencv/+cv/private/xxxx_.mexa64': Missing dependent shared libraries:
% 'libopencv_xfeatures2d.so.3.2' required by '/mexopencv/+cv/private/xxxx_.mexa64'
% read the Linux section of FAQ in this page : https://kyamagu.github.io/mexopencv/
% Simply, for the errors you may get using Kaze and KDtree, start matlab using:
% LD_PRELOAD=/home/tariq/Dev/opencv_build/lib/libopencv_xfeatures2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_features2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_flann.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_imgproc.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_core.so.3.2 matlab -nosplash -nodesktop
% You could also set an alias in your .bashrc:
% alias matlab='LD_PRELOAD=/home/tariq/Dev/opencv_build/lib/libopencv_xfeatures2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_features2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_flann.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_imgproc.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_core.so.3.2 matlab -nosplash -nodesktop'
MEXOPENCV_DIR = '/home/tariq/Dev/mexopencv' ; % '/path/to/mexopencv', for features extraction
VLFEAT_TOOLBOX_DIR = '/home/tariq/Dev/vlfeat/toolbox' ; % '/path/to/vlfeat/toolbox', for matching
MATLAB_CHOLMOD_DIR = '/home/tariq/Dev/SuiteSparse/CHOLMOD/MATLAB' ; % '/path/to/SuiteSparse/CHOLMOD/MATLAB', for spinv

% Simulation parameters :
FREQUENCY 		= 10 ; % frequency of acquisition for synchronisation
FIRST_IMAGE 	= 1001 ; % where to start reading the acquisition
LAST_IMAGE 		= 1100 ; % where to stop reading the acquisition
DOWN_SAMPLE 	= 3 ; % frames resampling frequency (next_frame = current_frame + DOWN_SAMPLE)
VIS_VERBOSE 	= 0 ; % show verbose during data acquisition ( 0-quiet, 1-text logging, 2-plot )

% Vision parameters
DETECTION_METHOD = 'KAZE' ; DETECTION_PARAMS = [0.001,10]; % KAZE: threshold,matching ratio between best two matches
%DETECTION_METHOD = 'FAST' ; DETECTION_PARAMS = [10,20,0,15]; % FAST: margins,threshold,non-max,LK-window
MIN_DISPARITY = 2.0 ; % minimum disparity to accept a correspondence
%RANSAC_ITR = 200 ; % RANSAC: number of iterations
%RANSAC_PIX_TOL = 0.5 ; % RANSAC: pixel tolerance
MIN_MUM_OF_EDGE_POINTS = 25 ; % minimum number of corners to include and edge
MIN_MUM_OF_EDGE_INLIERS = 15 ; % minimum number of outliers to compute a geometry

% optimiser parameters
NUM_OF_CAMS 	= 20 ; % PWGOPTIMISER: number of cams in each bundle
NUM_OF_VIEWS 	= 10 ; % PWGOPTIMISER: minimum number of views to accept a 3D point
IMAGE_NOISE 	= 1.0 ; % PWGOPTIMISER: image measurements noise
GATE_RATIO 		= 0.3 ; % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
MAX_NUM_ITR 	= 100 ; % PWGOPTIMISER: maximum number of iterations in case no convergence
OPT_VERBOSE 	= 1 ; % PWGOPTIMISER: show verbose during optimisation ( 0-quiet, 1-text logging, 2-plot )


