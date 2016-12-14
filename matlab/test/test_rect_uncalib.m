

% uncalibrated rectification
% Tariq Abuhashim - August 2014, iCub

function [H1, H2, im1_r, im2_r] = test_rect_uncalib(FM, cpts1, cpts2, im1, im2, im1g, im2g)

%
% notes:
%
% 1. Calculation of rectifying transformations H1 and H2 using both opencv and
%       matlab function works
% 2. Image rectification results were not good for disparity
%       estimation.
%

% number of points
npoints = size(cpts1, 1);
siz = size(im1g);

% collect points in cells points1 and points2
points1 = cell(npoints, 1);
points2 = points1;
for i = 1:npoints;
    points1{i} = cpts1(i, :);
    points2{i} = cpts2(i, :);
end

% plot epipolar lines before rectification
draw_epipolar(cpts1, cpts2, FM, im1, im2);

% H1 and H2 using opencv                             (this works very well)
% =================================================%
if 1
    [H1, H2] = cv.stereoRectifyUncalibrated(points1, points2, FM, siz);
    % check if H1 and H2 do the rectification job by plotting the
    % epipolar lines after rectification
    rhpts1 = H1*extend_points(cpts1'); rpts1 = flatten_points(rhpts1);
    rhpts2 = H2*extend_points(cpts2'); rpts2 = flatten_points(rhpts2);
    %F_rect = fundmatrix_nonlin([rpts1;rpts2], []);
    [F_rect, ~, ~] = test_Fmat(rpts1(1:2,:)', rpts2(1:2,:)', siz);
    disp('The rectified epipoles:');
    [e1, e2] = e1_e2_from_F(F_rect)
    draw_epipolar(rpts1', rpts2', F_rect, im1, im2);
end

% H1 and H2 using matlab                             (this works very well)
% =================================================%
if 0
    % T is the 3-by-3 transformation matrix that centres the image points
    origin = [siz(2); siz(1)]/2;
    T = [1 0 -origin(1); 0 -1 origin(2); 0 0 1];
    x1 = T*pextend(cpts1'); x2 = T*pextend(cpts2');
    [newim1, newim2, b, H1, H2] = rectify_uncalib(im1g, im2g, x1, x2, FM);
    
    % check if H1 and H2 do the rectification job by plotting the
    % epipolar lines after rectification
    rhpts1 = H1*extend_points(cpts1'); rpts1 = flatten_points(rhpts1);
    rhpts2 = H2*extend_points(cpts2'); rpts2 = flatten_points(rhpts2);
    %F_rect = fundmatrix_nonlin([rpts1;rpts2], []);
    [F_rect, ~, ~] = test_Fmat(rpts1(1:2,:)', rpts2(1:2,:)', siz);
    disp('The rectified epipoles:');
    [e1, e2] = e1_e2_from_F(F_rect)
    draw_epipolar(rpts1', rpts2', F_rect, newim1, newim2);
end

% warp the image using opencv          (rectification results are not good)
% =================================================%
if 0
    R1 = K_left\(H1*K_left);
    R2 = K_right\(H2*K_right);
    [mapx, mapy] = cv.initUndistortRectifyMap(K_left, [0 0 0 0 0], K_left, size(im1g'), 'R', R1);
    im1_r = cv.remap(im1, mapx, mapy);
    [mapx, mapy] = cv.initUndistortRectifyMap(K_right, [0 0 0 0 0], K_right, size(im2g'), 'R', R2);
    im2_r = cv.remap(im2, mapx, mapy);
    %dst1 = cv.warpPerspective(im1, H1);
    %dst2 = cv.warpPerspective(im2, H2);
    figure; imshow(im1_r); title('rect 1');
    figure; imshow(im2_r); title('rect 2');
    figure; imshow((im1_r+im2_r)/2);  title('sum of rects');
else
    im1_r = [];
    im2_r = [];
end

% using MATLAB to warp the image (this works better than opencv for the
% uncalibrated case, but the result images don't produce good disparity
% maps).                                   (disparity results are not good)
% =================================================%
if 1
    [t1, t2] = estimateUncalibratedRectification(FM, cpts1, cpts1, size(im2));
    tform1 = projective2d(t1);
    tform2 = projective2d(t2);
    % Rectify the images using projective transformations, tform1 and tform2.
    % Show a color composite of the rectified images demonstrating point
    % correspondences.
    I1Rect = imwarp(im1, tform1, 'OutputView', imref2d(size(im1)));
    I2Rect = imwarp(im2, tform2, 'OutputView', imref2d(size(im2)));
    % transform the points to visualize them together with the rectified images
    pts1Rect = transformPointsForward(tform1, cpts1);
    pts2Rect = transformPointsForward(tform2, cpts2);
    % show rectified images
    figure; showMatchedFeatures(I1Rect, I2Rect, pts1Rect, pts2Rect);
    legend('Inlier points in rectified I1', 'Inlier points in rectified I2');
    % Crop the overlapping area of the rectified images.
    % You can use red-cyan stereo glasses to see the 3D effect.
    [im_mosaic, im1_r, im2_r] = crop_and_transform(im1g, tform1, im2g, tform2);
    figure, imshow(im_mosaic);
    title('Rectified Stereo Images (Red - Left Image, Cyan - Right Image)');
    
    figure; imshow(im1_r); title('first rect image');
    figure; imshow(im2_r); title('second rect image');
else
    im1_r = [];
    im2_r = [];
end

% write images to files
% =================================================%
disp('writting rectified images to file ...')
croped_dist1=im1_r(20:end-20,20:end-20,:);
croped_dist2=im2_r(20:end-20,20:end-20,:);
imwrite(croped_dist1,...
'/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/img/qut_left.pgm');
imwrite(croped_dist2,...
'/home/tariq/Documents/MATLAB/koroibot/stereo/code/cpp/libelas/img/qut_right.pgm');