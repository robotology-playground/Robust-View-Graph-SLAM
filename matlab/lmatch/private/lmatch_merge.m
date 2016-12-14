% Merge the subsets of matches having a common segment in any image.
%
% l ... linesegs(:){1:K}, image line segments
% li ... double(K,:), set of matches that may contain common segments, indices to l
% ll ... lineseg(K,:), matched segments

function [MM,ll] = lmatch_merge(M,l,P,imsize)

K = size(M.li,1);

MM.li = [];

[lu{1:K}] = deal([]); % used segments
[ll{1:K}] = deal([]); % merged segments

while ~isempty(M.li)

  i = find_overlapping_matches(M.li,M.li(:,1));
  [lm,lui] = merge_overlapping_matches(l,M.li(:,i));
  
  j = size(MM.li,2) + 1;
  for k = 1:K
    if ~isempty(lm{k})
      ll{k} = [ll{k} lm{k}];
      lu{k} = [lu{k} lui{k}];
      MM.li(k,j) = length(ll{k});
    else
      MM.li(k,j) = 0;
    end
  end
  MM.C(j) = sum(M.C(i))/length(i);  % score of merged match is the mean of scores of original matches
  if i==1 % if no merging, keep the old L
    MM.X(:,j) = M.X(:,i);
    MM.Y(:,j) = M.Y(:,i);
  else  % if matches are really merged, we need to reconstruct the corresponding 3D line again
    %L = [sum([M.L(i).X],2) sum([M.L(i).Y],2)]/length(i); % initial estimate as the "mean" of original 3D lines
    k = any(M.li(:,i)>0,2);
    lmk = [lm{k}];
    L = line3d_from_lP_nonlin([lmk.u],[lmk.v],P(k),imsize(k),[],5,1e-8);
    L = lineseg3d_from_L(lmk,P(k),L);
    MM.X(:,j) = L(:,1);
    MM.Y(:,j) = L(:,2);
  end
  
  M.li(:,i) = [];
  M.C(i) = [];
  M.X(:,i) = [];
  M.Y(:,i) = [];
end

% add unused segments to ll
for k = 1:K
  ll{k} = [ll{k} l{k}(setminus(1:length(l{k}),lu{k}))];
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finds matches in li that have a common segment with n.
% li ... double(K,?)
% n ... double(K,1)
% i ... double(?), indices to li
function i = find_overlapping_matches(li,n)
i = [];
for k = 1:size(li,1)
  if n(k)~=0 % skip occluded/nondetected line segments
    i = [i find(n(k)==li(k,:))];
  end
end
i = unique(i);
return


% Merges matches li (indices to l) image-wise.
% lm ... lineseg(K), merged line segments; missing line segments are []
% lu ... cell{K} of double(?), used line segments, indices to l
function [lm,lu] = merge_overlapping_matches(l,li)
K = length(l);
[lu{1:K}] = deal([]);
for k = 1:K
  nk = unique(li(k,:));
  nk(nk==0) = []; % skip occluded/nondetected line segments
  if isempty(nk)
    lm{k} = [];
  else
    lm{k} = lineseg_merge(l{k}(nk));
    lu{k} = [lu{k} nk];
  end
end
return