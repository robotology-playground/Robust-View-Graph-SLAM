function options = check_data(DATA_DIR,robot)
%check_data(folder)
%
%	For all robots
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

% Robot name
% robot = "icub", "heicub" or "r1"
% Data related paths
% DATA_DIR = '/home/tariq/Documents/data/heicub/data_set1';

% Restore original <Default>Properties of root,
% load default PATH, run STARTUP.m:
% matlabrc;

% Robot related functions
addpath(strcat('./',robot)); % robot related functions

% Simulation parameters:
options = set_params(); % loads basic defaults
options = set_params(options,'folder',DATA_DIR); % where is the data?

[options, encoders, floatingbase] = set_images(options); % run data acquisition and synchronisation
													 % read inside notes to understand the 
													 % purpose of set_images
