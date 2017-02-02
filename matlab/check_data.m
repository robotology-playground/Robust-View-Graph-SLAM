function check_data(DATA_DIR,robot)
%check_data(folder)
%
%	For all robots
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

% Data related paths
%DATA_DIR='/home/tariq/Documents/data/heicub/data_set1';

% Restore original <Default>Properties of root,
% load default PATH, run STARTUP.m:
%matlabrc;

% Robot related functions
addpath(strcat('./',robot)); % robot related functions

% Simulation parameters:
options=set_params(); % loads basic defaults
% data collection parameters
options=set_params(options,'folder',		DATA_DIR); % where is the data?
options=set_params(options,'freq',			5	); 	% frequency of acquisition for synchronisation
													% Very large frequency means less accurate left to right synchronisation
													% Very small frequency means more frame drops
options=set_params(options,'verbose',		2	); 	% show verbose during data acquisition
													% 0 - no verbose, 1 - text logging, 2 - plot data
[options,encoders,floatingbase]=set_images(options); % run data acquisition and synchronisation
													 % read inside notes to understand the purpose of set_images

