function [H,B] = linesegs2H_metric(P1,P2,L)
% Compute homography H mapping neighborhood of segment u1 in image 1
% to neighborhood of segment u2 in image u2.
% There is 1-parameter family of homographies mapping l1 to l2 and consistent with cameras.
%
% The remaining parameter is computed using metric information, namely perpendicularity
% in the scene.
% H is induced by scene plane defined by scene line corresponding to the two line segments
% and by point at infinity being the normal of plane l1*P1 pulled back by the first segment.

L = vgg_contreps(L(:,1)*L(:,2)'-L(:,2)*L(:,1)');

X = hom(nhom(vgg_wedge(P1))+nhom(vgg_wedge(P2)));

A = X'*L;
B = [A(1:3) 0]*L;
H = P2*vgg_H_from_P_plane(B,P1);

return
