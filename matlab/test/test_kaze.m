function [kp1, d1, kp2, d2] = test_kaze(im1, im2, options)

% convert images
if size(im1, 3) > 1; 
    im1 = rgb2gray(im1); 
end;
if size(im2, 3) > 1; 
    im2 = rgb2gray(im2); 
end;

% extract features
[kp1, d1] = akaze(im1, 'dthreshold', options.kazethreshold, 'descriptor', 3);
fprintf([' ',num2str(size(kp1, 1)), ', ']);
[kp2, d2] = akaze(im2, 'dthreshold', options.kazethreshold, 'descriptor', 3);
fprintf([' ',num2str(size(kp2, 1)), '\n']);

if options.verbose > 1 % show features and images
    clf;
    subplot(1,2,1);imshow(im1);hold on;
    plot(kp1(:,1),kp1(:,2),'+');
    axis([1 size(im1,2) 1 size(im1,1)]);axis('ij');
    subplot(1,2,2);imshow(im2);hold on;
    plot(kp2(:,1),kp2(:,2),'+');
    axis([1 size(im2,2) 1 size(im2,1)]);axis('ij');
    drawnow;
end

% undistort points (no motion compensation)
%kp1 = remove_lens_distortion(kp1, options.kc1, options.K1);
%kp2 = remove_lens_distortion(kp2, options.kc2, options.K2);