
warning( 'off', 'MATLAB:mex:GccVersion' );

% use full paths, don't use ../ or ../../

vgslam = '/home/tariq/Documents/Robust-View-Graph-SLAM/'; 
eigen = '/home/tariq/Dev/eigen';
%suitesparse = '/home/tariq/Dev/suitesparse/'; % use dowloadable SuiteSparse
suitesparse = [vgslam 'SuiteSparse/']; % use supplied (old version) SuiteSparse

spinv = 1;
pwg = 1;
