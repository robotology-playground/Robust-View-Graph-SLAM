% Returns 1 iff two matches i(1) and i(2) are mergeable.
function b = lmatch_mergeable(i,M,l,P,I,D,opt)

% Construct merged line segments lm and test for merging residuals
for k = 1:length(l)
  if all(M.li(k,i)>0,2) & M.li(k,i(1))~=M.li(k,i(2))
    [b,lm(k)] = merge_test_2d(l{k}(M.li(k,i)),opt.MergeResidual);
    if ~b
      return
    end
  else
    maxn = max(M.li(k,i));
    if maxn>0
      lm(k) = l{k}(maxn);
    end
  end
end

% Leave only views that have segments in matches i
k = find(any(M.li(:,i)>0,2))';
lm = lm(k);
P = P(k);
I = I(k);
D = D(k,k);

% Reprojection test
[b,L] = lmatch_reprojection_test(lm,P,cellsizes(I),opt.ReprojResidual);
if ~b
  return
end

% Photommetric test
c = [];
for k = 1:length(P)
  [dummy,k1] = extreme_view(setminus(1:length(P),k),k,+1,D);
  c = [c min(lmatch_photoscore([lm(k).u lm(k).v],[lm(k1).u lm(k1).v],P{[k k1]},I{[k k1]},opt.NCCWindow,opt.Calibration>=3,L))];
end
b = all(c>opt.NCCThreshold(1));

return

    
% Returns 1 iff two linesegs are mergeable.
function [b,lm] = merge_test_2d(l,MergeResidual)
b = logical(0);
lm = l(1);

% test on scalar product of direction vectors
if sum(prod(normx([l.v]-[l.u])')) < .75
  return
end

% test on overlap
order12 = (l(2).v-l(2).u)'*(l(1).u-l(2).v) > 0; % 1 iff l(1) is before l(2)
order21 = (l(1).v-l(1).u)'*(l(2).u-l(1).v) > 0; % 1 iff l(2) is before l(1)
if ~( order12 | order21 )
  return
end

% fit residual test
s = [l.s];
lm = lineseg_merge(l);
lm.l = norml(cros(hom([lm.u lm.v]))');
ssd = vgg_vech(2*lm.l'*lm.l-diag(lm.l.^2))'*s;
if any(ssd./s(6,:) > MergeResidual^2)
  return
end
b = logical(1);

lm.s = sum(s,2)/length(s); % inv-covariance of merged segment; it is different from hom([lm.u lm.v])*hom([lm.u lm.v])'

return