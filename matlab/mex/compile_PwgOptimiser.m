%compile_PwgOptimiser( )
%
% A consule to compile and test mex files under construction
% This function uses the batch image optimisation class PwgOptimser and
% sparse inverse using Takahashi's method
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

config_install;

% COMPILE

if (spinv)
	spinv_install();	
end

if (pwg)
    pwg_install();
end
