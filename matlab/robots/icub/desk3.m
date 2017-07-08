% The following parameters were used in the FUSION2016 paper.
% Setting parameters for 20151207_desk3 dataset.
classdef desk3
	properties
		NAME					= 'desk3';
		% '/path/to/data'
		DATA_DIR 				= '/home/tariq/Documents/data/blue/20151207_desk3' ;
		% '/path/to/robot/calib_date'
		CALB_DIR 				= './icub/blue/calib_20151207'; 
		% Simulation parameters :
		FREQUENCY 				= 10; % frequency of acquisition for synchronisation
		FIRST_IMAGE 			= 51; % 51-where to start reading the acquisition
		LAST_IMAGE 				= 550; % 550-where to stop reading the acquisition
		DOWN_SAMPLE 			= 4; % frames resampling frequency (next_frame = current_frame + DOWN_SAMPLE)
		VIS_VERBOSE 			= 0; % show verbose during data acquisition ( 0-quiet, 1-text logging, 2-plot )
		% Vision parameters
		DETECTION_METHOD 		= 'FAST' ; 
		DETECTION_PARAMS 		= [10,20,0,15]; % FAST: margins,threshold,non-max,LK-window
		HORIZON 				= 50;
		MIN_DISPARITY 			= 10; % minimum disparity to accept a correspondence
		%RANSAC_ITR 			= 50 ; % RANSAC: number of iterations
		%RANSAC_PIX_TOL 		= 2 ; % RANSAC: pixel tolerance
		MIN_MUM_OF_EDGE_POINTS 	= 50; % minimum number of corners to include and edge
		MIN_MUM_OF_EDGE_INLIERS = 10; % minimum number of outliers to compute a geometry
		% optimiser parameters
		NUM_OF_CAMS 			= 20; % PWGOPTIMISER: number of cams in each bundle
		NUM_OF_VIEWS 			= 10; % PWGOPTIMISER: minimum number of views to accept a 3D point
		IMAGE_NOISE 			= 1.0; % PWGOPTIMISER: image measurements noise
		GATE_RATIO 				= 0.3; % PWGOPTIMISER: gate ratio between maximum and minimum acceptable inlier
		MAX_NUM_ITR 			= 100; % PWGOPTIMISER: maximum number of iterations in case no convergence
		OPT_VERBOSE 			= 1; % PWGOPTIMISER: show verbose during optimisation ( 0-quiet, 1-text logging, 2-plot )
	end
	methods
	end
end
