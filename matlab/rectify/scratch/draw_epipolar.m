

% Draws the epipolar lines on left and right images given matches and F
% Tariq Abuhashim - August 2014, iCub

function draw_epipolar(left_pts, right_pts, F, left_image, right_image)


% image dimensions
[~,n,~] = size(left_image);

% arrange data
left_x = left_pts(1,:); left_y = left_pts(2,:);
right_x = right_pts(1,:); right_y = right_pts(2,:);

% plot and show image on new figure
%figure;
subplot(2,1,1); imshow(left_image); axis image; hold on;
subplot(2,1,2); imshow(right_image); axis image; hold on;

% We know that left epipole is the 3rd column of V.
% We get V from svd of F. F=UDV'
[~, ~, FV] = svd(F);
left_epipole = FV(:,3);
left_epipole = left_epipole/left_epipole(3);
%right_epipole = FU(:,3);
%right_epipole = right_epipole/right_epipole(3);
    
% Start plotting:
for i=1:round(size(left_x,2)/10):size(left_x,2);
    
    % finding the epipolar line on the left image itself:
    % Hence using the left epipole and the given input point on left
    % image we plot the epipolar line on the left image
    left_epipolar_x = 1:n;
    left_epipolar_y = left_y(i) + (left_epipolar_x-left_x(i))*...
        (left_epipole(2)-left_y(i))/(left_epipole(1)-left_x(i));
    
    % plot on left image
    subplot(2,1,1);
    plot(left_epipolar_x, left_epipolar_y, 'r');
    plot(left_x(i), left_y(i), '*');
    
    % Getting the epipolar line on the RIGHT image (l_right = F x_left)
    % as a projection of the left point using the fundamental matrix
    left_P = [left_x(i); left_y(i); 1];
    right_P = F*left_P; % right epipolar line
    % Using the eqn of line: ax+by+c=0; y = (-c-ax)/b
    right_epipolar_x = 1:n;
    right_epipolar_y = (-right_P(3)-right_P(1)*right_epipolar_x)/right_P(2);
    
    % plot on right image
    subplot(2,1,2);
    plot(right_epipolar_x, right_epipolar_y, 'g');
    plot(right_x(i), right_y(i), '*');
    
end

drawnow;