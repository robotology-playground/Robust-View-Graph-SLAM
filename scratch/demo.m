clear all; close all; clc;

% parameter settings (for an example, please download
img_dir     = '/home/geiger/5_Data/kitti/2011_stereo/2010_03_09_drive_0019';
param.f     = 645.2;
param.cu    = 635.9;
param.cv    = 194.1;
param.base  = 0.571;
first_frame = 0;
last_frame  = 372;

% init visual odometry
StereoMex('init',param);

% read images
[im1, im2] = StereoMex('process',I1,I2);

% release visual odometry
%visualOdometryStereoMex('close');