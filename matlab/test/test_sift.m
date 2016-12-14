function [kp1, d1, kp2, d2] = test_sift(im1, im2, options)

% convert images
if size(im1,3) > 1; im1g = rgb2gray(im1); else im1g = im1; end;
if size(im2,3) > 1; im2g = rgb2gray(im2); else im2g = im2; end;
%if size(im1,3) > 1; im1g = im1(:,:,2); else im1g = im1; end;
%if size(im2,3) > 1; im2g = im2(:,:,2); else im2g = im2; end;

% % undistort images
% %%% Compute the new KK matx to fit as much data in the image (in order to
% %%% accomodate large distortions:
% %r2_extreme = (nx^2/(4*opt.fc1(1)^2) + ny^2/(4*opt.fc1(2)^2));
% dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
% fc_new = dist_amount * opt.fc1;
% KK_new = [fc_new(1) opt.alpha_c1*fc_new(1) opt.cc1(1);0 fc_new(2) opt.cc1(2) ; 0 0 1];
% im1g = rect(im1g,eye(3),opt.fc1,opt.cc1,opt.kc1,opt.alpha_c1,KK_new);
% im1g = im2uint8(im1g);
% %r2_extreme = (nx^2/(4*opt.fc2(1)^2) + ny^2/(4*opt.fc2(2)^2));
% dist_amount = 1; %(1+kc(1)*r2_extreme + kc(2)*r2_extreme^2);
% fc_new = dist_amount * opt.fc2;
% KK_new = [fc_new(1) opt.alpha_c2*fc_new(1) opt.cc2(1);0 fc_new(2) opt.cc2(2) ; 0 0 1];
% im2g = rect(im2g,eye(3),opt.fc2,opt.cc2,opt.kc2,opt.alpha_c2,KK_new);
% im2g = im2uint8(im2g);

% save undistorted image
%ima_name2 = [image_name '_rect.' format_image2];
%imwrite(I2,gray(256),ima_name2,format_image2);

Edge = options.siftthreshold(1);
Peak = options.siftthreshold(2);
% extract features
[kp1, d1] = vl_sift(im2single(im1g), 'EdgeThresh', Edge, 'PeakThresh', Peak);%,...
kp1 = double(kp1(1:2, :)');
d1 = double(d1);
fprintf([' ',num2str(size(kp1,1)),', ']);
[kp2, d2] = vl_sift(im2single(im2g), 'EdgeThresh', Edge, 'PeakThresh', Peak);
kp2 = double(kp2(1:2, :)');
d2 = double(d2);
fprintf([' ',num2str(size(kp2,1)),'\n']);

if options.verbose > 1 % show features and images
    clf;
    subplot(1,2,1);imshow(im1);hold on;
    plot(kp1(:,1),kp1(:,2),'+');
    axis([1 size(im1g,2) 1 size(im1g,1)]);axis('ij');
    subplot(1,2,2);imshow(im2);hold on;
    plot(kp2(:,1),kp2(:,2),'+');
    axis([1 size(im2g,2) 1 size(im2g,1)]);axis('ij');
    drawnow;
end