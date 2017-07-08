%
%	For HRP2@Toulouse.
%
% Tariq Abuhashim
% t.abuhashim@gmail.com
%
% Auckland, New Zealand, 2017

% Robot dataset :
fprintf(' - Using HRP2@Toulouse.\n') ;

param = circle1; % parameters for circle1 dataset

fprintf([' - Processing ',param.NAME,'.\n']) ;

% External functions :

% in order to compile mexopencv without installing opencv system wise:
% expert ${OPENCV_DIR}/lib/pkgconfig into PKG_CONFIG_PATH in your .bashrc 
% export PKG_CONFIG_PATH=${OPENCV_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH
% to make all, use: 	make all MATLABDIR=/usr/local/MATLAB/R2017a
% to make contrib, use: make contrib MATLABDIR=/usr/local/MATLAB/R2017a
% to make clean, use: 	make clean MATLABDIR=/usr/local/MATLAB/R2017a
% If you get an error: Invalid MEX-file '/mexopencv/+cv/private/xxxx_.mexa64': Missing dependent shared libraries:
% 'libopencv_xfeatures2d.so.3.2' required by '/mexopencv/+cv/private/xxxx_.mexa64'
% read the Linux section of FAQ in this page : https://kyamagu.github.io/mexopencv/
% Simply, for the errors you may get using Kaze and KDtree, start matlab using:
% LD_PRELOAD=/home/tariq/Dev/opencv_build/lib/libopencv_xfeatures2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_features2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_flann.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_imgproc.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_core.so.3.2 matlab -nosplash -nodesktop
% You could also set an alias in your .bashrc:
% alias matlab='LD_PRELOAD=/home/tariq/Dev/opencv_build/lib/libopencv_xfeatures2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_features2d.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_flann.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_imgproc.so.3.2:/home/tariq/Dev/opencv_build/lib/libopencv_core.so.3.2 matlab -nosplash -nodesktop'

MEXOPENCV_DIR 		= '/home/tariq/Dev/mexopencv' ; % '/path/to/mexopencv', for features extraction
VLFEAT_TOOLBOX_DIR 	= '/home/tariq/Dev/vlfeat/toolbox' ; % '/path/to/vlfeat/toolbox', for matching
MATLAB_CHOLMOD_DIR 	= '/home/tariq/Dev/SuiteSparse/CHOLMOD/MATLAB' ; % '/path/to/SuiteSparse/CHOLMOD/MATLAB', for spinv

% end of file
