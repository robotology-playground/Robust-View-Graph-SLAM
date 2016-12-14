function [im1, im2] = get_image(e1, e2, options)

if mod(e1, 2);
    idx = (e1+1)/2;
    im1 = imread(strcat(options.img_folder{1},options.cam_left.image{idx}));
else
    idx = e1/2;
    im1 = imread(strcat(options.img_folder{2},options.cam_right.image{idx})); 
end;

if mod(e2, 2);
    idx = (e2+1)/2;
    im2 = imread(strcat(options.img_folder{1},options.cam_left.image{idx}));
else
    idx = e2/2;
    im2 = imread(strcat(options.img_folder{2},options.cam_right.image{idx})); 
end;

% bayer decoding?
if size(im1, 3) == 1
    fprintf('Bayer decoding: ');
    im1 = demosaic(im1, 'grbg');
    im2 = demosaic(im2, 'grbg');
end