% Reprojection test on K line segments in K images
%
% l ... lineseg(K), line segments in K images
% L ... double(4,2), 3D line previously reconstructed from l(1:K-1) (can be omitted)
% P ... cell(K)
% imsize ... double(2,K)
% reproj_residual ... reprojection residual threshold (in pixels)
% b ... logical, result of the test
%
function [b,L] = lmatch_reprojection_test(l,P,imsize,reproj_residual,L)

K = length(P);
Pj = vertcat(P{:});

% Covariance matrices of line segments
s = [l.s];
for k = 1:K
  ss{k} = vgg_vech(s(:,k));
end

% If residuals after linear estimation are very large, reject
Llin = vgg_line3d_from_lP_lin(ss,P,imsize);
ll = cros(reshape(Pj*Llin(:,1),[3 K]),reshape(Pj*Llin(:,2),[3 K]))';
ssd = sum(s.*vgg_vech_swap(norml(ll)'));
if max(ssd./s(6,:)) > reproj_residual(1)^2
  b = logical(0);
  L = zeros(1,6);
  return
end

% Otherwise, do non-linear estimation
if nargin < 5
  L = Llin;
end
if 1 % faster
  L = line3d_from_lP_nonlin([l.u],[l.v],P,imsize,L,5,1e-8);
else % slower
  L = vgg_line3d_from_lP_nonlin(ss,P,imsize,L,[],'niter_term',4,'loglambda_init',-15,'verbose',0);
end
ll = cros(reshape(Pj*L(:,1),[3 K]),reshape(Pj*L(:,2),[3 K]))';
ssd = sum(s.*vgg_vech_swap(norml(ll)'));
b = max(ssd./s(6,:)) < reproj_residual(2)^2; % reprojection residual threshold

return