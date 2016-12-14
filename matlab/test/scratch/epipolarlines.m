function [leftx,lefty,rightx,righty] = plot_epipolar_lines(left_x, left_y, FM, im1, im2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function takes as argument the (x,y) coordinate of the LEFT
% image and plots its corresponding epipolar line on the RIGHT and  
% the left image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard coding the fundamental matrix FM

% FM=[  0.0001   -0.0437   13.3897
%       0.0429    0.0033  -21.4157
%      -12.1830   22.0779 -999.3629];

% close all;
% left_image  = double(rgb2gray(imread('rodin0020.jpg')));
% right_image = double(rgb2gray(imread('rodin0021.jpg')));

[m,n,~] = size(im1);


% Start plotting:

% Getting the epipolar line on the RIGHT image:
left_P = [left_x; left_y; 1];
right_P = FM*left_P;
rightx = 1:n;
% Using the eqn of line: ax+by+c=0; y = (-c-ax)/b
righty = (-right_P(3)-right_P(1)*rightx)/right_P(2);
% plot
subplot(1,2,1); imshow(im1); axis image; hold on;
plot(left_epipolar_x, left_epipolar_y, 'r');

% Now finding the other epipolar line on the left image itself:
% We know that left epipole is the 3rd column of V.
% We get V from svd of F. F=UDV'
[~, ~, FV] = svd(FM);
left_epipole = FV(:,3);
left_epipole = left_epipole/left_epipole(3);
% Hence using the left epipole and the given input point on left
% image we plot the epipolar line on the left image
leftx = 1:n;
lefty = left_y + (leftx-left_x)*(left_epipole(2)-left_y)/(left_epipole(1)-left_x);
% plot
subplot(1,2,2); imshow(im2); axis image; hold on;
plot(left_epipolar_x, left_epipolar_y, 'g');
