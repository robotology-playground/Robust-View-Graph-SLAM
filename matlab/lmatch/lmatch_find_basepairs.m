% Ki = lmatch_find_basepairs(D,N)  Finds disjoint basepairs "most evenly spread" across the view set.
% D ... double(#views,#views), pairwise view distance (see lmatch_generate)
% N ... required number of basepairs

function Ki = lmatch_find_basepairs(D)

K = size(D,1);

% find two closest views
for k = 1:K
  [d(k),i(k)] = extreme_view(setminus(1:K,k),k,+1,D);
end
k = argmin(d);
Ki = [k; i(k)];

while 1
  % k := farthest view from existing set Ki
  d = [];
  for k = setminus(1:K,Ki(:))
    [d(k),i(k)] = extreme_view(Ki(:),k,+1,D);
  end
  if isempty(d)
    break
  end
  k = argmax(d);

  % find nearest view to k
  nKi = setminus(1:K,[Ki(:)' k]);
  if isempty(nKi)
    break
  end
  [dummy,i] = extreme_view(nKi,k,+1,D);
  Ki = [Ki [k;i]];
end

return
