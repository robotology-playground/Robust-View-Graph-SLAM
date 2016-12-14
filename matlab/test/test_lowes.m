function [kp1,d1,kp2,d2]=test_lowes(im1,im2,options)

if size(im1,3) > 1; im1g = rgb2gray(im1); else im1g = im1; end;
if size(im2,3) > 1; im2g = rgb2gray(im2); else im2g = im2; end;
% extract sift
[desc, kpts] = mod_sift_lowe(im1g);
kp1 = kpts(:, [2 1]);
d1 = desc';
[desc, kpts] = mod_sift_lowe(im2g);
kp2 = kpts(:, [2 1]);
d2 = desc';

if options.verbose % show features
    clf;
    subplot(1,2,1); imshow(im1); hold on;
    plot(kp1(:,1), kp1(:,2), '+');
    subplot(1,2,2); imshow(im2); hold on;
    plot(kp2(:,1), kp2(:,2), '+');
    drawnow;
end