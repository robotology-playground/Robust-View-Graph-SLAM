function img = test_undistort_image(img, fc, alpha_c, cc, kc)

if size(img, 3) > 1;
    img = rgb2gray(img);
 end

        
%%% Compute the new KK matrix to fit as much data in the image (in order to
%%% accomodate large distortions:
[ny, nx] = size(img);
r2_extreme = (nx^2/(4*fc(1)^2) + ny^2/(4*fc(2)^2));
dist_amount = 1; % removes empty pixels (image is smaller)
%dist_amount = (dist_amount + kc(1)*r2_extreme + kc(2)*r2_extreme^2); % keeps empty pixels
fc_new = dist_amount * fc;
KK_new = [fc_new(1) alpha_c*fc_new(1) cc(1); 0 fc_new(2) cc(2) ; 0 0 1];
img = rect(img, eye(3), fc, cc, kc, alpha_c, KK_new);
img = im2uint8(img);