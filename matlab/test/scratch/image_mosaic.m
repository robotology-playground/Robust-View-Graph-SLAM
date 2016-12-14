% 
% Produces a mosaic of images im1 and im2 gives their projective
% transformation H
%
%
% Tariq Abuhashum
% started: 21 August 2014
%
% iCub - Koroibot
%

function mosaic=image_mosaic(im1,im2,H)

% mosaic
box2 = [1  size(im2,2) size(im2,2)  1 ;
    1  1           size(im2,1)  size(im2,1) ;
    1  1           1            1 ] ;
box2_ = H\box2 ;
box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;

[u,v] = meshgrid(ur,vr) ;
im1_ = vl_imwbackward(im2double(im1),u,v) ;

z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
im2_ = vl_imwbackward(im2double(im2),u_,v_) ;

mass = ~isnan(im1_) + ~isnan(im2_) ;
im1_(isnan(im1_)) = 0 ;
im2_(isnan(im2_)) = 0 ;
mosaic = (im1_ + im2_) ./ mass ;

figure ; clf ;
imagesc(mosaic) ; axis image off ;
title('Mosaic') ;

if nargout == 0;
    clear mosaic;
end