%
% plot_features.m
%
% Tariq Abuhashim
% started: 19 August 2014
%
% iCub - Koroibot
%

function test_plot_features(im1, im2, kpts1, kpts2, mpts1, mpts2, cpts1, cpts2)

% dimensions
[w, h, c] = size(im1);

% start
figure(1);

% raw features
subplot(2,2,1); hold off;
imshow(im1); hold on;
plot(kpts1(:,1), kpts1(:,2),'g+');
title([num2str(size(kpts1, 1)), ' features']);
subplot(2,2,2); hold off;
imshow(im2); hold on;
plot(kpts2(:,1), kpts2(:,2),'r+');
title([num2str(size(kpts1, 1)), ' features']);
drawnow;

% matched features
subplot(2,2,3);
showMatchedFeatures(im1, im2, mpts1', mpts2');
legend('matched points 1','matched points 2');
title([num2str(size(mpts1, 2)),' matches']);
drawnow;

% inliers (after mapsac, rotation, and visibility)
subplot(2,2,4);
showMatchedFeatures(im1, im2, cpts1', cpts2');
legend('inliers 1','inliers 2');
title([num2str(size(cpts1, 2)),' inliers']);

drawnow;