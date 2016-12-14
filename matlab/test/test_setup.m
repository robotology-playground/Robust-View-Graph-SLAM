

% setup folders, parameters and read images
% Tariq Abuhashim - August 2014, iCub

function [im1, im2, im1g, im2g] = test_setup( )


% code
addpath('/home/tariq/Documents/MATLAB/koroibot/sba/code');
addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/matlab/common');
addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/matlab/rectify');
% features
addpath('/home/tariq/Dev/vlfeat-0.9.16/toolbox'); vl_setup;
addpath('/home/tariq/Dev/akaze-master/mex');
% external resources
addpath('/home/tariq/Documents/MATLAB/koroibot/external/torrsam');
addpath('/home/tariq/Documents/MATLAB/koroibot/external/mexopencv');
% disparity
addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/matlab');
addpath('/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/src');


% read images

% image folder and base names
% folder{1} = '/home/tariq/Documents/MATLAB/koroibot/stereo/data/UQ_St_Lucia/';
% folder{2} = '/home/tariq/Documents/MATLAB/koroibot/stereo/data/UQ_St_Lucia/';
% im1_name = 'cam0_image01520.png';
% im2_name = 'cam1_image01520.png';
folder{1} = '/home/tariq/Documents/MATLAB/koroibot/stereo/data/data_test_1/img/left/';
folder{2} = '/home/tariq/Documents/MATLAB/koroibot/stereo/data/data_test_1/img/right/';
im1_name = '00001999.ppm';
im2_name = '00001999.ppm';
im1 = imread([folder{1}, im1_name]);
im2 = imread([folder{2}, im2_name]);

% make grayscale
[~, ~, depth] = size(im1);
if depth > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
if depth > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end