
%
% The Epipolar line search - plane/parallax decomposition
%
% section 3 of paper: Motion - Stereo Integration for Depth Estimation
% Christopher Strecha and Luc Van Gool
% ECCV 2002
%
% Chapter 13 of book: Multiple View Geometry in Computer Vision
% Richard Hartley and Andrew Zisserman
% second edition 2003
%
% Tariq S. Abuhashim
%

function d = epipolargeometry(F, t, R, x1, x2, CameraParams)

K  = [CameraParams.fu 0               CameraParams.u0;
      0               CameraParams.fv CameraParams.v0;
      0               0               1];

d = -100:.1:100; % depth ambiguity

hfig_epipolar = figure;

subplot(1,2,1);
axis([0 CameraParams.umax 0 CameraParams.vmax]);
grid on;
axis equal;
hold on;

subplot(1,2,2);
axis([0 CameraParams.umax 0 CameraParams.vmax]);
grid on;
axis equal;
hold on;

for itr = 1:500:size(x1, 2);
    
    % first image
    xl = [x1(1,itr); x1(2,itr); 1];
    figure(hfig_epipolar);
    subplot(1,2,1);
    axis([0 CameraParams.umax 0 CameraParams.vmax]);
    plot(xl(1), xl(2), 'O');
    
    % second image (dense stereo)
    xr = zeros(2, length(d));
    for i=1:length(d)
        H_inf = K*R'/K; % homography of plane at infinity
        line = H_inf*xl - d(i)*K*R'*t;
        xr(1,i) = line(1)/line(3);
        xr(2,i) = line(2)/line(3);
    end
    figure(hfig_epipolar);
    subplot(1,2,2);
    axis([0 CameraParams.umax 0 CameraParams.vmax]);
    plot(xr(1,:), xr(2,:));
    
    xr_h = (K*R'/K)*xl; % correspondence assuming planner scene
    xr(1) = xr_h(1)/xr_h(3);
    xr(2) = xr_h(2)/xr_h(3);
    plot(xr(1), xr(2),'O');
    
    imagesize = [CameraParams.umax CameraParams.vmax];
    plot(x2(1,itr), x2(2,itr),'rO');
    [leftx,lefty,rightx, righty] = epipolarlines(x2(1,itr), x2(2,itr), F, imagesize);
    plot(leftx,lefty, 'r');
    
%     % second image (Motion parallax)
%     x = xl(1) - K(1,3);
%     y = xl(2) - K(2,3);
%     Qr = [x*y/K(1,1)        -x*x/K(1,1)-K(1,1) y;
%           y*y/K(2,2)+K(2,2) -x*y/K(2,2)       -x];
%     Qt = [K(1,1) 0      -x;
%           0      K(2,2) -y];
%     xr = zeros(2, length(d));
%     for i=1:length(d)
%         xr(:,i) = xl(1:2) + Qr*Cb2n2euler(R)' + d(i)*Qt*t;
%     end
%     figure(hfig_epipolar);
%     subplot(1,2,2);
%     axis([0 CameraParams.umax 0 CameraParams.vmax]);
%     plot(xr(1,:), xr(2,:),'r');
%     epipolar = xl(1:2) + Qr*Cb2n2euler(R)';
%     plot(epipolar(1), epipolar(2),'gO');
%     drawnow;
end

%close(hfig_epipolar);

%
% END
%