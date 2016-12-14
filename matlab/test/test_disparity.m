

% disparity estimation
% Tariq Abuhashim - August 2014, iCub

function [D1, D2] = test_disparity(dst1, dst2, opt)

% read images
%images_folder='/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas';
%im1=imread([images_folder,'/img/qut_left.pgm']);
%im2=imread([images_folder,'/img/qut_right.pgm']);
im1 = dst1;im2 = dst2;

% convert images
if size(im1,3) > 1; im1g = rgb2gray(im1); else im1g = im1; end;
if size(im2,3) > 1; im2g = rgb2gray(im2); else im2g = im2; end;

% disparity estimation parameters
my_params ( );

% estimate the disparity
[D1, D2] = elasMex(im1g', im2g', param);
%[D1,D2] = mex_desc(im1g',im2g',param);