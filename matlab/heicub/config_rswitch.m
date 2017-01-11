% function config_rswitch()
%function config_rswitch()
%
% Loads residuals switch configurations
%
%	For iCub@heidelberg.
%
% Tariq Abuhashim - 2016
% t.abuhashim@gmail.com
%
% iCub - Koroibot


sigma_r = 1/options.K1(1,1); % pixels
R = eye(2)*sigma_r*sigma_r;

terr=0.001; % norm of translation error
aerr=0.0002; % norm of rotation error
delta=[terr, aerr]; % delta errors
gmin=0.30; % minimum residuals gate probability
gmax=0.95; % maximum residuals gate probability
maxitr=options.maxitr; % maximum number of optimisation iterations

%GATE_PROB_RESID=gmin+rand(1,maxitr)*(gmax-gmin);

%GATE_PROB_RESID = rand(1, maxitr);
%GATE_PROB_RESID = ones(1,maxitr) * .5;
%GATE_PROB_RESID = min(GATE_PROB_RESID, .95); % max. prob. limit
%GATE_PROB_RESID = max(GATE_PROB_RESID, .05); % min. prob. limit
%GATE_PROB_RESID(1) = max(GATE_PROB_RESID(1),.5); % min. first gate prob. limit

GATE_PROB_RESID=[0.05 0.9 0.7]; % residuals gating probability of inliers

TRUST_INIT_LINKS=false; % if set TRUE, the set Ct = C(sw==1);
