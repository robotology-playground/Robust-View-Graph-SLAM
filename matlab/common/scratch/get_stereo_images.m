
function [cam_left, cam_right] = get_stereo_images(folder, base, s_1, s_2, s_3)

% [cam_left, cam_right] = get_stereo_images(folder, base, s_1, s_2, s_3)
% folder : folder where images reside (can be separate for left and right images)
% base : the base of image names (examples: 'cam0_image*.png' or  '*.ppm')
% s_1 : first image number, (default = 1)
% s_2 : last image number, (default = number of images in the folder)
% s_3 : image step (jumps), (default = 1)
%
% get_stereo_images.m
%
% Tariq Abuhashim
% started: 6 August 2014
%
% iCub - Koroibot

% check if left and right camera images are in separate folders
if iscell(folder) && size(folder,2) > 1;
    left_folder = folder{1};
    right_folder = folder{2};
else
    left_folder = folder;
    right_folder = folder;
end
if iscell(base) && size(base, 2) > 1;
    left_base = base{1};
    right_base = base{2};
elseif iscell(base);
    left_base = base{1};
    right_base = base{1};
else
    left_base = base;
    right_base = base;
end

% get image names
filename_1 = dir(strcat(left_folder, left_base));
filename_2 = dir(strcat(right_folder, right_base));
n_1 = size(filename_1, 1); % total number of images
n_2 = size(filename_2, 1);
disp(['there are ',num2str(n_1), ' and ', ...
    num2str(n_2), ' images in camera one and two folders.']);

% check image read settings
if nargin < 3; s_1 = 1; end; % start
if nargin < 4; s_2 = n_1; end; % stop
if nargin < 5; s_3 = 1; end; % step
disp(['running images ',num2str(s_1),' to ',...
    num2str(s_2),' in steps of ',num2str(s_3),'.']);
if s_2 > n_1; 
    error('img_stop is larger than the number of images in the folder'); end

% camera one (left) names
k = 0;
for i = s_1:s_3:s_2;
    k = k+1;
    cam_left.image{k} = filename_1(i).name;
end

% camera two (right) names
k = 0;
for i = s_1:s_3:s_2;
    k = k+1;
    cam_right.image{k} = filename_2(i).name;
end