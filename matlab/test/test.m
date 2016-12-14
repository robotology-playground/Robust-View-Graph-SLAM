

% suplement code can be found at:
% Andrew Zisserman : http://www.robots.ox.ac.uk/~vgg/hzbook/code/
% Peter Kovesi : http://www.csse.uwa.edu.au/~pk/research/matlabfns/
% Jean-Yves Bouguet : http://www.vision.caltech.edu/bouguetj/calib_doc/
%
% Tariq Abuhashim - August 2014, iCub

clc; clear all; close all;

[im1, im2, im1g, im2g] = test_setup( );

% image processing
%[im1, im2, im1g, im2g] = test_image_processing(im1, im2);

% 2d points
[kpts1, desc1, kpts2, desc2] = test_features(im1g, im2g, 'kaze');

% 2d matching
[mpts1, mpts2, matches] = test_matching(kpts1, kpts2, desc1, desc2);

% fundamental matrix                                         (robustness?)
[FM, cpts1, cpts2] = test_Fmat(mpts1, mpts2, size(im1g));

% essential matrix
[P1, P2, R_est, t_est, U_est, vis] = test_Emat(cpts1, cpts2);
R_est
t_est

% check features
figure;
test_plot_features(im1, im2, kpts1, kpts2, mpts1, mpts2, cpts1, cpts2);

% check epipolar lines
figure;
draw_epipolar(cpts1, cpts2, FM, im1, im2);

% plot 3d points
figure;
plot3(U_est(1,:), U_est(2,:), U_est(3,:), '.');
axis equal; grid on;

return;


% uncalibratedrectification        (for points, not for images, no scale)
[H1, H2, im1_r, im2_r] = test_rect_uncalib(FM, cpts1, cpts2, im1, im2, im1g, im2g);


% calibrated rectification            (not for points, for images, scale?)
[R1, R2, T1, T2, dst1, dst2] = test_rect_calib(FM, cpts1, cpts2, im1, im2, im1g, im2g);


% disparity estimation                               (robustness?, scale?)
D = test_disparity(dst1, dst2);


% computer depth
[d1, d2] = test_depth( );
% plots
test_plot_disparity;


% jacobean
xf = [1 2 3; 4 5 6; 7 8 9]; 
xc = [2 3 4 .1 .1 .1]';
[H, e] = test_jacobian(xf, xc);