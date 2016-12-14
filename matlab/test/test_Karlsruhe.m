% test Karlsruhe dataset

close all;
clc;
clear all;

% folders
addpath('/home/tabuhashim/Dev/akaze/mex');
addpath('/home/tabuhashim/Dev/vlfeat-0.9.19/toolbox'); vl_setup;
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/torrsam');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/mexopencv');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/common');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/rectify');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/data/Karlsruhe/2011_09_26');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot//stereo/code/cpp/libelas/matlab');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot//stereo/code/cpp/libelas/src');

% images
base_dir = '/home/tabuhashim/Documents/MATLAB/koroibot/stereo/data/Karlsruhe';
folder{1} = [base_dir, '/2011_09_26/2011_09_26_drive_0052_extract/image_00/data/'];
folder{2} = [base_dir, '/2011_09_26/2011_09_26_drive_0052_extract/image_01/data/'];
% folder{1} = [base_dir, '/2011_09_26/2011_09_26_drive_0052_sync/image_02/data/'];
% folder{2} = [base_dir, '/2011_09_26/2011_09_26_drive_0052_extract/image_03/data/'];

folder{1} = [base_dir, '/2009_09_08/2009_09_08_drive_0010/left/'];
folder{2} = [base_dir, '/2009_09_08/2009_09_08_drive_0010/right/'];
base{1} = 'I1_*.png';
base{2} = 'I2_*.png';
[cam_left, cam_right] = get_stereo_images(folder, base, 1, 200, 1);

% calibration
%calibration;
roi = [1.800000e+01 1.361000e+03 6.000000e+01 4.500000e+02];
K1 = [9.032949e+02 0.000000e+00 6.639935e+02 0.000000e+00 9.079042e+02 2.452070e+02 0.000000e+00 0.000000e+00 1.000000e+00];
D1 = [-3.778799e-01 1.824904e-01 1.390637e-03 4.659340e-05 -4.730213e-02];
R1 = [9.998650e-01 -1.241631e-02 -1.076237e-02 1.239131e-02 9.999204e-01 -2.387250e-03 1.079115e-02 2.253568e-03 9.999392e-01];
K2 = [9.050234e+02 0.000000e+00 6.846818e+02 0.000000e+00 9.102276e+02 2.457892e+02 0.000000e+00 0.000000e+00 1.000000e+00];
D2 = [-3.810763e-01 1.851087e-01 1.802898e-03 -1.702626e-04 -4.624721e-02];
R2 = [9.996878e-01 -1.176894e-02 -2.203912e-02 1.182005e-02 9.999277e-01 2.190397e-03 2.201175e-02 -2.450217e-03 9.997547e-01];
R = [9.999369e-01 -5.437215e-04 1.122319e-02 5.966196e-04 9.999887e-01 -4.710484e-03 -1.122050e-02 4.716883e-03 9.999259e-01];
T = [-5.724638e-01 6.739395e-03 1.262054e-02];
P1 = [6.790081e+02 0.000000e+00 6.778034e+02 0.000000e+00 0.000000e+00 6.790081e+02 2.465724e+02 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];
P1_roi = [6.790081e+02 0.000000e+00 6.598034e+02 0.000000e+00 0.000000e+00 6.790081e+02 1.865724e+02 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];
P2 = [6.790081e+02 0.000000e+00 6.778034e+02 -3.887481e+02 0.000000e+00 6.790081e+02 2.465724e+02 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];
P2_roi = [6.790081e+02 0.000000e+00 6.598034e+02 -3.887481e+02 0.000000e+00 6.790081e+02 1.865724e+02 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00];

K1 = reshape(K1, 3, 3)'; %D1 = D00;
K2 = reshape(K2, 3, 3)'; %D2 = D01;
Rc = reshape(R, 3, 3)'; %Rc = Rc';
%Rc = reshape(R_01, 3, 3)'\reshape(R_00, 3, 3)'; % better results
tc = T'; %T_01' - T_00';
P1_rect = reshape(P1, 4, 3)';
P2_rect = reshape(P2, 4, 3)';

% stereo reconstruction
for i = 100;%:size(cam_right.image, 2);
    
    % rectification
    im1 = imread([folder{1},cam_left.image{i}]);
    im2 = imread([folder{2},cam_right.image{i}]);
    if size(im1, 3)>1; im1g = rgb2gray(im1); else im1g = im1; end;
    if size(im2, 3)>1; im2g = rgb2gray(im2); else im2g = im2; end;
    S = cv.stereoRectify(K1, D1, K2, D2, size(im1g'), Rc, tc);
    % S.R1– 3x3 rectification transform (rotation matrix) for the first camera.
    % S.R2– 3x3 rectification transform (rotation matrix) for the second camera.
    % S.P1– 3x4 projection matrix in the new (rectified) coordinate systems for the first camera.
    % S.P2– 3x4 projection matrix in the new (rectified) coordinate systems for the second camera.
    % S.Q – 4x4 disparity-to-depth mapping matrix (reprojectImageTo3D).
    
%     e1 = tc/norm(tc);
%     dz = [0 0 1]'; e2 = cross(e1, dz); e2 = e2/norm(e2);
%     e3 = cross(e1, e2);
%     Rr = [e1 e2 e3]'; R1 = Rr; R2 = Rc*Rr;
    
    [mapx, mapy] = cv.initUndistortRectifyMap(K1, D1, P1_rect, size(im1g'));%,'R', S.R1);
    dst1 = cv.remap(im1, mapx, mapy);
    [mapx, mapy] = cv.initUndistortRectifyMap(K2, D2, P2_rect, size(im1g'));%,'R', S.R2);
    dst2 = cv.remap(im2, mapx, mapy);
    %box = [max([S.roi1(1:2); S.roi2(1:2)]) min([S.roi1(3:4); S.roi2(3:4)])];
    %dst1 = dst1(box(2):box(4),box(1):box(3),:);
    %dst2 = dst2(box(2):box(4),box(1):box(3),:);
%     figure(1); clf; imshow(dst1); title('rect 1');
%     figure(2); clf; imshow(dst2); title('rect 2');
    figure(3); clf; imshow((dst1+dst2)/2); 
    title('sum of rects');
    
    % check the epipolar lines of rectified images
    [k1, d1, k2, d2] = test_features(dst1, dst2, 'kaze');
    [m1, m2, mat] = test_matching(k1, k2, d1, d2, 3);
    [F, cpts1, cpts2] = test_Fmat(m1, m2, size(im1g), .005, 1000, .001);
    null(F)
    figure(4); clf; draw_epipolar(cpts1, cpts2, F, dst1, dst2);
    
%     % disparity
%     [D1, D2] = test_disparity(dst1, dst2);
%     figure(5); clf; imagesc(D1');
%     figure(5); clf; imagesc(D2');
%     
%     % 3d from disparity (not actual depth map)
%     test_plot_disparity;
    
    %pause;
    
end