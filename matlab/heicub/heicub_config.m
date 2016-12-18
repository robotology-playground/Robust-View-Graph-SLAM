function [options,encoders,floatingbase]=heicub_config()
% [options,encoders,floatingbase]=heicub_config()
%
% Tariq Abuhashim
% started: September 2014
%
% iCub - Koroibot

DATA_DIR='/home/tariq/Documents/data/heicub/data_set1';
CALB_DIR='/home/tariq/Documents/data/heicub/calib_20160913/img/stereo';
SAVE_DIR=strcat(DATA_DIR,'/run_',datestr(now,'yyyymmdd'));
 
options=set_params(); % loads basic defaults
options=set_params(options,'folder',		DATA_DIR);
options=set_params(options,'save',			SAVE_DIR);
options=set_params(options,'calib',			CALB_DIR);
options=set_params(options,'freq',			10);
options=set_params(options,'first_image',	50); % 50
options=set_params(options,'last_image',	80); % 800
options=set_params(options,'steps',			2	);
options=set_params(options,'verbose',		0	);
 
[options,encoders,floatingbase]=set_images(options); % get the images using options
 
options=set_params(options,'vision'); % loads the vision defaults
options=set_params(options,'fastthreshold',	20.0);
options=set_params(options,'mindisp',		2.0	);
options=set_params(options,'ransac',		200	);
options=set_params(options,'RANSAC_pixtol',	0.5	);
options=set_params(options,'mincorrnr',		50	);
options=set_params(options,'mininlnr',		25	);
options=set_params(options,'gridmargin',	5.0	);
options=set_params(options,'gridhorizon',	5.0	);
 
options=set_params(options,'optimiser'); % loads the pwg_optimiser defaults
options=set_params(options,'ncams',			5	);
options=set_params(options,'nkeys',			1	);
options=set_params(options,'nview',			10	);
options=set_params(options,'sigma_r',		1.0	);
options=set_params(options,'gateratio',		0.2	);
options=set_params(options,'verbose',		1	);
options=set_params(options,'maxitr',		100	);
