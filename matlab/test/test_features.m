

% extract features
% Tariq Abuhashim - August 2014, iCub

function [kpts1, desc1, kpts2, desc2] = test_features(im1, im2, method)

% convert images
if size(im1,3) > 1; im1g = rgb2gray(im1); else im1g = im1; end;
if size(im2,3) > 1; im2g = rgb2gray(im2); else im2g = im2; end;

% Kaze features
if strcmp(method, 'kaze');
    
    fprintf('kaze ... '); t_start = tic;
    [kpts1, desc1] = akaze(im1g, 'dthreshold', .0001, 'descriptor', 3);
    [kpts2, desc2] = akaze(im2g, 'dthreshold', .0001, 'descriptor', 3);
    
elseif strcmp(method, 'sift');

    fprintf('sift ... '); t_start = tic;
    [kpts1, desc1] = vl_sift(im2single(im1g)); 
    kpts1 = kpts1(1:2, :)';
    %desc1 = single(desc1);
    [kpts2, desc2] = vl_sift(im2single(im2g)); 
    kpts2 = kpts2(1:2, :)';
    %desc2 = single(desc2);
    
elseif strcmp(method, 'dsift');

    fprintf('dsift ... '); t_start = tic;
    [kpts1, desc1] = vl_dsift(im2single(im1g), 'size', 20); 
    kpts1 = kpts1(1:2, :)';
    [kpts2, desc2] = vl_dsift(im2single(im2g), 'size', 20); 
    kpts2 = kpts2(1:2, :)';
    
elseif strcmp(method, 'cvsift');

    fprintf('cvsift ... '); t_start = tic;
    [feat, desc1] = cv.SIFT(im1g); 
    kpts1 = zeros(size(feat, 2), 2); desc1 = desc1';
    for i = 1:size(feat, 2);
        kpts1(i, :) = feat(i).pt;
    end
    [feat, desc2] = cv.SIFT(im2g); 
    kpts2 = zeros(size(feat, 2), 2); desc2 = desc2';
    for i = 1:size(feat, 2);
        kpts2(i, :) = feat(i).pt;
    end 
    
elseif strcmp(method, 'cov');
    
    fprintf('cov ... '); t_start = tic;
    [kpts1, desc1] = vl_covdet(im2single(im1g), 'EstimateOrientation', true);
    kpts1 = kpts1(1:2, :)';
    [kpts2, desc2] = vl_covdet(im2single(im2g), 'EstimateOrientation', true);
    kpts2 = kpts2(1:2, :)';
    
end

% print text
fprintf([num2str(size(kpts1,1)),' left, ',num2str(size(kpts2,1)),' right, ']);
t_end=toc(t_start); fprintf([num2str(t_end),'s\n']);