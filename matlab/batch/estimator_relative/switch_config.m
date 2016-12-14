% % Configurations
% USE_VISION = 0;
% TRUST_INIT_LINKS = 0;
% options.sigma_r = 2; % 2 pixels
% options.gateratio = 0.7;
% options.gateinnov = chi_square_bound(.9,2);
% GATE_PROB_RESID = [0.9,0.9,0.7,0.5];
% SIGMA_0 = 1e+06;
% SIGMA_t = 1e+06;
% SIGMA_a = 1e+02;
% SIGMA_p = 1e-04;

% Configurations
USE_VISION = 1; % 0-kinematics initialisation, 1-vision initialisation
TRUST_INIT_LINKS = 0;
STEREO = 0; % calibrated stereo is assumed between pair {[1,2], [3,4], ...}
%GATE_PROB_RESID = [.5,.5,.5,.5,.5,.5,.5,.5,.5,.5];

%options.sigma_r = 0.5; % 0.5 pixels
%options.gateratio = 0.2;
%options.gateinnov = chi_square_bound(.99,2); %.99
%iSIGMA2_0 = 1/(1e-4)^2; % 0.01mm
%iSIGMA2_t = 1/(2.5e-3)^2; % 5.00mm
%iSIGMA2_a = 1/(3*pi/180)^2;% large error
%iSIGMA2_p = 1/(1e+5)^2; % 10000mm
%options.verbose = 2;