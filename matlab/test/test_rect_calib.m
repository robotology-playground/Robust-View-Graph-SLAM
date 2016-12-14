

% calibrated rectification
% Tariq Abuhashim - August 2014, iCub

function [R1, R2, T1, T2, dst1, dst2] = test_rect_calib(FM, cpts1, cpts2, im1, im2, im1g, im2g)

%
% notes:
%
% 1. Need to implement the calculation of rectifying transformations H1,H2
%       where xhat1=H1*x1 and xhat2=H2*x2, and epipolar lines are parallel
% 2. Image rectification results were good for disparity estimation.
%

% number of points
npoints = size(cpts1, 1);

% camera parameters
[K_left, K_right, d_left, d_right, Rc, tc] = test_param( );

% collect points in cells points1 and points2
points1 = cell(npoints, 1);
points2 = points1;
for i = 1:npoints;
    points1{i} = cpts1(i, :);
    points2{i} = cpts2(i, :);
end

% plot epipolar lines before rectification
draw_epipolar(cpts1, cpts2, FM, im1, im2);

% calibrated rectification using opencv
S = cv.stereoRectify(K_left, d_left, K_right, d_right, size(im1g), Rc, tc);

% outputs
R1 = S.R1; R2 = S.R2; 
T1 = S.P1; T2 = S.P2;

% plot epipolar lines before rectification ?????????? (this is to be added)

% warp the image using opencv                        (this works very well)
% =================================================%
[mapx, mapy] = cv.initUndistortRectifyMap(K_left, d_left, T1(:,1:3), size(im1g'),'R', R1);
dst1 = cv.remap(im1, mapx, mapy);
[mapx, mapy] = cv.initUndistortRectifyMap(K_right, d_right, T2(:,1:3), size(im1g'),'R', R2); 
dst2 = cv.remap(im2, mapx, mapy);
figure; imshow(dst1); title('rect 1');
figure; imshow(dst2); title('rect 2');
figure; imshow((dst1+dst2)/2); title('sum of rects');

% write images to files
% =================================================%
disp('writting rectified images to file ...')
croped_dist1=dst1(20:end-20,20:end-20,:);
croped_dist2=dst2(20:end-20,20:end-20,:);
imwrite(croped_dist1,...
'/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/img/qut_left.pgm');
imwrite(croped_dist2,...
'/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/img/qut_right.pgm');