% this file tests the functionality of P1_P2_from_F

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
[P1_est, P2_est, R_est, t_est, U_est, vis, e_est] = test_Emat(cpts1, cpts2);

%% tests
clc
F = F_from_P1_P2(P1_est, P2_est);
E = E_from_R_t(R_est, t_est);
[P1, P2] = P1_P2_from_F(E);

P2_est
P2