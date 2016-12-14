% 
% Performs image mosaicing of images im1 and im2 using sift features and
% optimises the projective transformation between the two images
%
% Tariq Abuhashum
% started: 21 August 2014
%
% iCub - Koroibot
%

function mosaic = sift_mosaic(im1, im2)

% load folder
vlfeat_folder='/home/tariq/Dev/vlfeat-0.9.16/toolbox';
addpath(vlfeat_folder); vl_setup;

% read image
if nargin == 0
    im1 = imread(fullfile(vl_root, 'data', 'river1.jpg')) ;
    im2 = imread(fullfile(vl_root, 'data', 'river2.jpg')) ;
end

% make single
im1 = im2single(im1) ;
im2 = im2single(im2) ;

% make grayscale
if size(im1,3) > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
if size(im2,3) > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end

% features matching
[f1,d1] = vl_sift(im1g) ;
[f2,d2] = vl_sift(im2g) ;
[matches, scores] = vl_ubcmatch(d1,d2) ;
numMatches = size(matches,2) ;

% computer homography
[H,inliers]=homography_ransac(f1(1:2,matches(1,:)),f2(1:2,matches(2,:)));

% show matches
show_matches(im1,im2,f1(1:2,:),f2(1:2,:),matches,inliers);

% mosaic
mosaic=image_mosaic(im1,im2,H);

end