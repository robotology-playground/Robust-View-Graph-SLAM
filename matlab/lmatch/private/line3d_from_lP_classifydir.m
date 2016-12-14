% [L,d] = line3d_from_lP_classifydir(s,P,imsize,L0,V,ReprojResid)  Reconstructs 3D line, simultaneously classifying it into a direction class.
%
% s, P, imsize ... see vgg_line3d_from_lP_nonlin
% V ... cell(N), list of direction classes:
%   - If V{d} is double(4,1), L passes through homog point V{d}.
%   - If V{d} is double(4,2), L passes through line spanned by homog points V{d}(:,1) and V{d}(:,2).
%   - If V{d}==[], L is unconstrained.
% Lines with 2, 3, and 4 DOF are tried in turn. If either has reprojection residual smaller then ReprojResid,
% it is accepted.

function [L,d] = line3d_from_lP_classifydir(s,P,imsize,L0,V,ReprojResid)

if nargin < 6
  ReprojResid = 1.5;
end
if nargin < 5
  ReprojResid = inf;
  V = {[]};
end
if nargin < 4
  L0 = [];
end

K = length(P);

% for d = 1:length(V)
%   L = vgg_line3d_from_lP_nonlin(s,P,imsize,[],V{d},'niter_term',10,'drmsrel_term',1e-9);
%   for k = 1:K
%     l(k,:) = norml(cros(P{k}*L)');
%     ssd(k) = l(k,:)*s{k}*l(k,:)';
%     n(k) = s{k}(3,3);
%   end
%   rms(d) = sqrt(sum(ssd)/sum(n));
% end
% disp(rms)
% return

p = cellsizes(V);
p = p(2,:);
q{1} = find(p==1);
q{2} = find(p==2);
q{3} = find(p==0);

for f = 1:3
  if isempty(q{f})
    continue
  end
  rms = [];
  L = [];
  for d = q{f}
    L(:,:,d) = vgg_line3d_from_lP_nonlin(s,P,imsize,L0,V{d});
    for k = 1:K
      l(k,:) = norml(cros(P{k}*L(:,:,d))');
      ssd(k) = l(k,:)*s{k}*l(k,:)';
      n(k) = s{k}(3,3);
    end
    rms(d) = sqrt(sum(ssd)/sum(n));
  end

  [r,i] = min(rms(q{f}));
  if r < ReprojResid
    d = q{f}(i);
    L = L(:,:,d);
    break
  end
end

if r >= ReprojResid
  d = 0;
end

return