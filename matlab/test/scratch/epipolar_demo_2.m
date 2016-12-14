
function epipolar_demo()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function shows a demo of epipolar geometry.
%
% The user should click a point in the left image. The corresponding
% epipolar line for the clicked point is shown on right image.(The
% corresponding epipolar line on the left image is also shown.)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading the correspondence points in the two image:

load left_image_points;
load right_image_points;

[a b] = size(left_image_points);

% finding A matrix:

for i=1:a
  
    x1 = left_image_points(i,1);
    y1 = left_image_points(i,2);
    x2 = right_image_points(i,1);
    y2 = right_image_points(i,2);
    A(i,:) = [x1*x2 y1*x2 x2 x1*y2 y1*y2 y2 x1 y1 1];
    
end

% SVD of A:
[U D V] = svd(A);

% Finding Fundamental Matrix F:
f = V(:,9);
F = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];

% Modify F:
[FU FD FV]= svd (F);
FDnew = FD;
FDnew(3,3) = 0;

FM = FU*FDnew*FV';

save FM;

% Plotting epipolar line:

close all;
left_image = double(rgb2gray(imread('rodin0020.jpg')));
right_image = double(rgb2gray(imread('rodin0021.jpg')));

[m n] = size(left_image);

figure,imagesc(left_image); colormap(gray); title('Click a point on this Left Image');axis image;
figure,imagesc(right_image); colormap(gray); title('Corresponding Epipolar Line in this Right Image');axis image;

list =['r' 'b' 'g' 'y' 'm' 'k' 'w' 'c'];

for i=0:7
    
    % Clicking a point on the left image:
    figure(1);    
    [left_x left_y] = ginput(1);
    hold on;
    plot(left_x,left_y,'r*');

    % Finding the epipolar line on the right image:
    left_P = [left_x; left_y; 1];

    right_P = FM*left_P;

    right_epipolar_x=1:2*m;
    % Using the eqn of line: ax+by+c=0; y = (-c-ax)/b
    right_epipolar_y=(-right_P(3)-right_P(1)*right_epipolar_x)/right_P(2);
    figure(2);
    hold on;
    plot(right_epipolar_x,right_epipolar_y,list(mod(i,8)+1));

    % Now finding the other epipolar line
    % We know that left epipole is the 3rd column of V.
    % We get V from svd of F. F=UDV'
    left_epipole = FV(:,3);
    left_epipole = left_epipole/left_epipole(3);
    
    left_epipolar_x = 1:2*m;
    left_epipolar_y = left_y + (left_epipolar_x-left_x)*(left_epipole(2)-left_y)/(left_epipole(1)-left_x);
    figure(1);
    hold on;
    plot(left_epipolar_x,left_epipolar_y,list(mod(i,8)+1));

end
