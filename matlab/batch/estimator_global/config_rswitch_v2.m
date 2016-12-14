%%% Configuration file

% for compute_switches algorithm
% GATE_PROB_INNOV = 0.99;
% GATE_PROB_RESID = [0.05 0.3 0.9 0.7]; % for version 2, 0.05 is designed to force cooperation with Ct
% GATE_PROB_RESID = [0.9 0.9 0.7]; % for version 3
% K_RESID = 0.3;

% % not bad with run_20151015
% USE_MAX_SPANNING_TREE = false;
% NUM_INITIAL_PATH = 5;
% REFINE_INITIAL_PATH = true; % only used if ~USE_MAX_SPANNING_TREE && NUM_INITIAL_PATH > 1
% TRUST_INIT_LINKS = false;
% INIT_WITH_HIGH_WEIGHT_LINKS = 0;
%
% GATE_PROB_INNOV = 0.8;
% GATE_PROB_RESID = [0.3 0.9 0.7];
% K_RESID = 0.3;

USE_MAX_SPANNING_TREE = 1;
NUM_INITIAL_PATH = 3;
REFINE_INITIAL_PATH = 1; % only used if ~USE_MAX_SPANNING_TREE && NUM_INITIAL_PATH > 1
TRUST_INIT_LINKS = 0;
INIT_WITH_HIGH_WEIGHT_LINKS = 0;

GATE_PROB_INNOV = .99;
%GATE_PROB_RESID = [0.5, 0.3, 0.05];
GATE_PROB_RESID = [0.3, 0.9, 0.7, 0.5];
K_RESID = 0.3; %0.9

% % when using kinematics
% SIGMA_t = 1/1000;
% SIGMA_a = .2*pi/180;

% % when being pushed
% SIGMA_t = 5/1000;
% SIGMA_a = 1*pi/180;

% when using RGBD
SIGMA_t = 5/1000; % in meters
SIGMA_a = 1*pi/180; % in radians

VERBOSE = 1;

EDGE_WEIGHT = 1;
% wighting method
% 1 - loop statistics (RANSAC)
% 2 - number of tracks
% 3 - trace of information
% 4 - trace of covariance
