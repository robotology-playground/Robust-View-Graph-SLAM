
clc; clear all; close all;

folder = '/home/tabuhashim/Documents/MATLAB/koroibot/stereo/data/ikea/black/head_1';
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/common');
addpath('/home/tabuhashim/Documents/MATLAB/koroibot/stereo/code/matlab/patch');
opt.FOLDER=folder;
opt.IMAGE_BASE='*.ppm';
opt.baseline=1;
opt.freq=30;
opt.STEPS=1;
opt.FIRST_IMAGE=670;
opt.LAST_IMAGE=730;
opt.save='/home/tabuhashim/Documents/MATLAB/koroibot/stereo/data/ikea/black/results/head_1/run_11';
opt = patch_opt( opt );
set_folders;
opt = set_params( opt );
[cam_left, cam_right, eyes, neck, waist, img_folder] = set_images(opt);


for i = 1;%:size(cam_left.image,2)
    
    % read images
    im1g = rgb2gray(im2double(imread([img_folder{1}, cam_left.image{i}])));
    im2g = rgb2gray(im2double(imread([img_folder{2}, cam_right.image{i}])));
    
%     % undistort images
%     im1g = test_undistort_image(im1g,opt.imgsize(2),opt.imgsize(1), ...
%     opt.fc1,opt.alpha_c1,opt.cc1,opt.kc1);
%     im2g = test_undistort_image(im2g,opt.imgsize(2),opt.imgsize(1), ...
%     opt.fc2,opt.alpha_c2,opt.cc2,opt.kc2);
    
    % mosaic and estimate homography
    [mosaic, H, im1_, im2_] = sift_mosaic(im1g, im2g);
    
%     % rectification
%     S = cv.stereoRectify(opt.K1, opt.kc1, opt.K2, opt.kc2, size(im1g'), eye(3), [.068; 0; 0]);
%     [mapx, mapy] = cv.initUndistortRectifyMap(opt.K1, opt.kc1, S.P1, size(im1g'),'R', S.R1\S.R1);
%     dst1 = cv.remap(im1g, mapx, mapy);
%     [mapx, mapy] = cv.initUndistortRectifyMap(opt.K2, opt.kc2, S.P2, size(im2g'),'R', S.R2\S.R1);
%     dst2 = cv.remap(im2g, mapx, mapy);
    
    % mosaic the rectified images
    mass=~isnan(dst1)+~isnan(dst2) ;
    dst1(isnan(dst1))=0;
    dst2(isnan(dst2))=0;
    imshow((dst1 + dst2)./mass);

    pause;
    
end