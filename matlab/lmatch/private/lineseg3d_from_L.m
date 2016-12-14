% S = lineseg3d_from_L(s,P,L)  Reconstructs single 3D line segment.
%
% s ... lineseg(K)
% P ... cell(K) of double(3,4)
% L ... double(4,2), 3D line
% X ... double(4,2), homog. end points of the 3D segment

function X = lineseg3d_from_L(s,P,L)

K = length(P);

% project segment's end points on reprojected straight lines
Pj = vertcat(P{:});
l = vgg_wedge(reshape(Pj*L(:,1),[3 K]),reshape(Pj*L(:,2),[3 K])); % reprojected image lines
for k = 1:K
  Q{k} = subtx(l(k,:));
  x{k} = Q{k}*subinv(Q{k})*hom([s(k).u s(k).v]);
end

% re-project all end points into image 1
x1 = x{1};
for k = 2:K
  x1 = [x1 vgg_contreps(l(1,:))*vgg_F_from_P(P{[k 1]})*x{k}];
end

if 0
Lpm = L(:,1)*L(:,2)' - L(:,2)*L(:,1)';
lnu = [l(:,2) -l(:,1)];
lnv = [l(:,2) -l(:,1)];
X = [];
for k = 1:K
  lnu(k,3) = -lnu(k,1:2)*x{k}(1:2,1)/x{k}(3,1);
  lnv(k,3) = -lnv(k,1:2)*x{k}(1:2,2)/x{k}(3,2);
  ln = [lnu(k,:);lnv(k,:)];
  X = [X (ln*P{k}*Lpm)'];
end
end

% x := most distant points of x1
t = nhom(subinv(Q{1})*x1);
t = [min(t) max(t)];
t0 = nhom(subinv(Q{1})*hom([s(1).u s(1).v]));
if prod([t; t0]*[1;-1]) < 0 % set correct order of end points to preserve segment orientation
  t = t([2 1]);
end
x = Q{1}*hom(t);

% reconstruct 3D end points
A = l(2,:)*P{2};
X = normx([vgg_contreps(P{1}'*vgg_contreps(x(:,1))*P{1})*A',...
           vgg_contreps(P{1}'*vgg_contreps(x(:,2))*P{1})*A']);

return