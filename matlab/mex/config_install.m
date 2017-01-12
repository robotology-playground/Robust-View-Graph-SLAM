%config_install( )
%
% A consule to compile and test mex files under construction
% This function loads compile configurations
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Koroibot, iCub Facility, Istituto Italiano di Tecnologia
% Genova, Italy, 2016

% Disable mex compiler compatability warning.
warning( 'off', 'MATLAB:mex:GccVersion' );

% Use full paths, don't use ../ or ../../
vgslam = '/home/tariq/Documents/Robust-View-Graph-SLAM/'; 
eigen = '/home/tariq/Dev/eigen';
%eigen = '/usr/local/include/eigen*'; % eigen or eigen3
%suitesparse = '/home/tariq/Dev/suitesparse/'; % use dowloadable SuiteSparse
suitesparse = [vgslam 'SuiteSparse/']; % use supplied (old version) SuiteSparse

spinv = 1;
pwg = 1;
