function kpts = remove_lens_distortion(kpts, k, K)

%%% the distortion model is:

% xd and yd are the "distorted" version of the need coordinates x and y
% xd = x.*(1+k1*r2 + k2*r2.^2) + 2*p1.*x.*y + p2*(r2 + 2*x.^2);
% yd = y.*(1+k1*r2 + k2*r2.^2) + 2*p2.*x.*y + p1*(r2 + 2*y.^2);
rev = 0;
if size(kpts, 2) > 2;
    kpts = kpts';
    rev = 1;
end

k1 = k(1);
k2 = k(2);
p1 = k(3);
p2 = k(4);
k3 = k(5);
if nargin > 2
    fx = K(1,1);
    fy = K(2,2);
    cx = K(1,3);
    cy = K(2,3);
    xd = (kpts(:,1) - cx)/fx;
    yd = (kpts(:,2) - cy)/fy;
else % when using normalised coordinates input/output
    xd = kpts(:,1);
    yd = kpts(:,2);
end
x = xd;
y = yd;
%for i = 1:size(xd, 1);
%   [x(i), y(i)] = UndistPointInv(xd(i), yd(i), k1, k2, k3, p1, p2);
%end
[x, y] = UndistPointInv(xd, yd, k1, k2, k3, p1, p2);
% remove normalisation
if nargin > 2
    x = fx*x + cx;
    y = fy*y + cy;
end
% outputs
kpts(:,1:2) = [x y];
if rev;
    kpts = kpts';
end

%%% undistort points using MATLAB symbolic toolbox
% syms xd yd k1 k2 p1 p2 r2 x y
% equ1 = x.*(1+k1*r2+k2*r2.^2)+2*p1.*x.*y+p2*(r2+2*x.^2)==xd;
% equ2 = y.*(1+k1*r2+k2*r2.^2)+2*p2.*x.*y+p1*(r2+2*y.^2)==yd;
% %eqns = [equ1, equ2];
% S = solve(equ1, equ2, x, y)


%r2 = xd.^2 + yd.^2;
% x = -((xd - p2.*r2).*(yd - 2.*yd.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + ...
%     (k1.*r2)/2 + (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - ...
%     2.*p2.^2.*r2 + 2.*xd.*p2 + 1/4).^(1/2) - p1.*r2 - k1.*p1.*r2.^2 - k2.*p1.*r2.^3 + ...
%     k2.*r2.^2.*yd + 2.*p1.*r2.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + ...
%     (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + ...
%     2.*xd.*p2 + 1/4).^(1/2) + k1.*r2.*yd))./(4.*(yd - p1.*r2).*(- r2.*p1.^2 + ...
%     yd.*p1 - r2.*p2.^2 + xd.*p2));
%
% y = -(yd - 2.*yd.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + ...
%     (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + ...
%     2.*xd.*p2 + 1/4).^(1/2) - p1.*r2 - k1.*p1.*r2.^2 - k2.*p1.*r2.^3 + k2.*r2.^2.*yd + ...
%     2.*p1.*r2.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + (k2.^2.*r2.^4)/4 + ...
%     (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + 2.*xd.*p2 + 1/4).^(1/2) + ...
%     k1.*r2.*yd)./(4.*(- r2.*p1.^2 + yd.*p1 - r2.*p2.^2 + xd.*p2));


%r2 = xd.^2 + yd.^2;
% x = -((xd - p2.*r2).*(yd + 2.*yd.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + ...
%     (k1.*r2)/2 + (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - ...
%     2.*p2.^2.*r2 + 2.*xd.*p2 + 1/4).^(1/2) - p1.*r2 - k1.*p1.*r2.^2 - k2.*p1.*r2.^3 + ...
%     k2.*r2.^2.*yd - 2.*p1.*r2.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + ...
%     (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + ...
%     2.*xd.*p2 + 1/4).^(1/2) + k1.*r2.*yd))./(4.*(yd - p1.*r2).*(- r2.*p1.^2 + ...
%     yd.*p1 - r2.*p2.^2 + xd.*p2));
% 
% y = -(yd + 2.*yd.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + ...
%     (k2.^2.*r2.^4)/4 + (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + ...
%     2.*xd.*p2 + 1/4).^(1/2) - p1.*r2 - k1.*p1.*r2.^2 - k2.*p1.*r2.^3 + k2.*r2.^2.*yd - ...
%     2.*p1.*r2.*((k1.^2.*r2.^2)/4 + (k1.*k2.*r2.^3)/2 + (k1.*r2)/2 + (k2.^2.*r2.^4)/4 + ...
%     (k2.*r2.^2)/2 - 2.*p1.^2.*r2 + 2.*yd.*p1 - 2.*p2.^2.*r2 + 2.*xd.*p2 + 1/4).^(1/2) + ...
%     k1.*r2.*yd)./(4.*(- r2.*p1.^2 + yd.*p1 - r2.*p2.^2 + xd.*p2));