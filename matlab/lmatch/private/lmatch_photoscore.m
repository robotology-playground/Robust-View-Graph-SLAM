% c = lmatch_photoscore(u1,u2,P1,P2,I1,I2,CorrParams,metric [,L])
%
% Returns similarity score of two line segments for individual pixels.
%
% u1, u2 ... double (2,2), line segment in image 1 resp. 2.
% P1, P2 ... cameras, oriented
% I1, I2 ... double (:,:), images
% CorrParams ... correlation window parameters, [wr_orthog wr_parall pixdist both_sides]
% c ... vector, x-correlations of line points
%
% If metric==0, oriented projective calibration of cameras is assumed.
% If metric==1, metric calibration is assumed.
%
% L ... double(4,2), 3D infinite line - stabilizes geometry if L was reconstructed from 3+ views
%
% Line segments are clipped by epipolars inside the function.
% Line segments must be consistently oriented, their directions matter.

function c = lmatch_photoscore(u1,u2,P1,P2,I1,I2,CorrParams,metric,L)

if nargin < 9
  L = [];
end

if isempty(L) % reconstruct L from image lines
  l1 = cros(hom(u1))';
  l2 = cros(hom(u2))';
  L = null([l1*P1; l2*P2]);
else
  l1 = cros(P1*L)';
  l2 = cros(P2*L)'; % L projected to images
end

% X := segment in image 2 projected onto L
n = [-l2(2) l2(1)];
X2 = ( [n -n*u2(:,1); n -n*u2(:,2)]*P2*(L(:,1)*L(:,2)'-L(:,2)*L(:,1)') )';

% Clip segment in image 1 by reprojected segment from image 2
Q = subtx(hom(u1));
iQ = subinv(Q);
t1 = sort(nhom(iQ*hom(u1)));
t2 = sort(nhom(iQ*P1*X2));
t = boxmeet(t1,t2);
if any(isnan(t))
  c = 0;
  return
end
v1 = nhom(Q*hom(t));
if [1 -1]*u1'*v1*[1;-1] > 0
  u1 = v1;
else
  u1 = v1(:,[2 1]);
end
%t = nhom(subinv(Q)*[hom(u1) P1*X2]);
%t = sort(t);
%t = t([3 2]); % This is done because we know that t(1)>t(2), which follows from det(R)>0 where subtx(x)*R=x. Thus, orientation of the segment is preserved.
%u1 = nhom(Q*hom(t));

% homography mapping of neighborhoods of the line segments
if metric
  H = linesegs2H_metric(P1,P2,L);
else
  H = linesegs2H_oriproj(nhom(P1*X2),u2,P1,P2);
end

c = scoreH(u1,H,I1,I2,CorrParams);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function c = scoreH(u1,H,I1,I2,CorrParams)

% correlation window displacements
wr = CorrParams(1:2); % wr = halfsides of correlation window, in respectively orthogonal and parallel direction to the line segment
i = ones(2*wr(2)+1,1)*[-wr(1):wr(1)];
j = (ones(2*wr(1)+1,1)*[-wr(2):wr(2)])';
kk{1} = 1:2*wr(1)+1;
kk{2} = 1:wr(1)+1;
kk{3} = wr(1)+1:2*wr(1)+1;

% line segment pixels
pixdist = CorrParams(3); % distance of neighboring corr. window centers on the line
p1 = uv2pixels(u1,pixdist);

% rotation to line direction
l1 = cros(hom(u1))';
d1 = normx(l1(1:2)');
R1 = [d1 [0 -1; 1 0]*d1];

% Loop for correlating respectively middle, left and right neighborhood of the line segment.
% Of these 3 scores, the best is selected in each segment pixel.
c = -ones(1,size(p1,2));
switch CorrParams(4)
 case 0, K = 1;
 case 1, K = 1:3;
end
for k = K
  ki = i(:,kk{k});
  kj = j(:,kk{k});
  w = R1*[ki(:)'; kj(:)'];
  w(3,:) = 1;
  c = max(c,vgg_ncc_2im_H(I1,I2,H,w,p1));
end

return


% x = uv2pixels(uv,step)  Given end points of a line segment, returns pixels from u to v with step.
function x = uv2pixels(u,step)
du = u*[-1;1];
p = du/sqrt(du'*du);
t = p'*du;
t = 0 : step*sign(t) : t;
x = p*t + u(:,ones(1,size(t,2)));
return