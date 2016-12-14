function [rect_im1,rect_im2,H]=rect_right_to_left(im1,im2,x1,x2)

% this function is my version of the function sift_mosaic.m
% it projects the second image into the first image by estimating 
% the mapping homography with RANSAC.

% make grayscale
if size(im1,3) > 1, im1g = rgb2gray(im1) ; else im1g = im1 ; end
if size(im2,3) > 1, im2g = rgb2gray(im2) ; else im2g = im2 ; end

% check if points are in homog-coords
if size(x1,1)>size(x1,2);
    x1=x1'; x2=x2';
    
end
numMatches = size(x1,2) ;
X1 = x1(1:2,:); X1(3,:) = 1;
X2 = x2(1:2,:); X2(3,:) = 1;

% RANSAC with homography model
clear H score ok ;
for t = 1:100
    % estimate homograpyh
    subset = vl_colsubset(1:numMatches, 4) ;
    A = [] ;
    for i = subset
        A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
    end
    [U,S,V] = svd(A) ;
    H{t} = reshape(V(:,9),3,3) ;
    
    % score homography
    X2_ = H{t} * X1 ;
    du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
    dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
    ok{t} = (du.*du + dv.*dv) < 6*6 ;
    score(t) = sum(ok{t}) ;
end

[score, best] = max(score) ;
H = H{best} ;
ok = ok{best} ;

% Optional refinement
    function err = residual(H)
        u = H(1) * X1(1,ok) + H(4) * X1(2,ok) + H(7) ;
        v = H(2) * X1(1,ok) + H(5) * X1(2,ok) + H(8) ;
        d = H(3) * X1(1,ok) + H(6) * X1(2,ok) + 1 ;
        du = X2(1,ok) - u ./ d ;
        dv = X2(2,ok) - v ./ d ;
        err = sum(du.*du + dv.*dv) ;
    end

if exist('fminsearch') == 2
    H = H / H(3,3) ;
    opts = optimset('Display', 'none', 'TolFun', 1e-8, 'TolX', 1e-8) ;
    H(1:8) = fminsearch(@residual, H(1:8)', opts) ;
else
    warning('Refinement disabled as fminsearch was not found.') ;
end

% Show matches
dh1 = max(size(im2,1)-size(im1,1),0) ;
dh2 = max(size(im1,1)-size(im2,1),0) ;

figure(1) ; clf ;
subplot(2,1,1) ;
imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
o = size(im1,2) ;
line([x1(1,:);x2(1,:)+o], [x1(2,:);x2(2,:)]) ;
title(sprintf('%d tentative matches', numMatches)) ;
axis image off ;

subplot(2,1,2) ;
imagesc([padarray(im1,dh1,'post') padarray(im2,dh2,'post')]) ;
o = size(im1,2) ;
line([x1(1,ok);x2(1,ok)+o], [x1(2,ok);x2(2,ok)]) ;
title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
    sum(ok), ...
    100*sum(ok)/numMatches, ...
    numMatches)) ;
axis image off ;

drawnow ;

% --------------------------------------------------------------------
%                                                               Mosaic
% --------------------------------------------------------------------

box2 = [1  size(im2,2) size(im2,2)  1 ;
    1  1           size(im2,1)  size(im2,1) ;
    1  1           1            1 ] ;
box2_ = inv(H) * box2 ;
box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;

[u,v] = meshgrid(ur,vr) ;
rect_im1 = vl_imwbackward(im2double(im1),u,v) ;

z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
rect_im2 = vl_imwbackward(im2double(im2),u_,v_) ;

mass = ~isnan(rect_im1) + ~isnan(rect_im2) ;
im1_(isnan(rect_im1)) = 0 ;
im2_(isnan(rect_im2)) = 0 ;
mosaic = (rect_im1 + rect_im2) ./ mass ;

figure(2) ; clf ;
imagesc(mosaic) ; 
axis image off ;
title('Mosaic') ;

end
