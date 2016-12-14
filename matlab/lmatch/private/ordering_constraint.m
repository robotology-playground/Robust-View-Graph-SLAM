% N = ordering_constraint(li0,li,l,P)  Subset of li incompatible with li0 due to ordering constraint.
%
% li0 ... double(K,1), query match, indices to l
% li ... double(K,?), matching match, indices to l
% l ... cell(K) of linesegs, image line segments
% P ... cell{K}, cameras, correctly oriented
% i ... indices of incompatible matches

function i = ordering_constraint(i0,ii,M,l,P,alg)

i = [];
for k1 = 1:length(P)
  ii1 = ii( M.li(k1,ii)>0 & M.li(k1,i0)~=M.li(k1,ii) ); % matches with existing segment and having no segment in common with i0
  for k2 = k1+1:length(P)
    k = [k1 k2];
    if any(M.li(k,i0)<=0) % query match must have segments in both views k1,k2
      continue
    end
    ii2 = ii1( M.li(k2,ii1)>0 & M.li(k2,i0)~=M.li(k2,ii1) ); % matches with existing segment and having no segment in common with i0
    %ii2 = ii( all(M.li(k,ii)>0) & M.li(k1,i0)~=M.li(k1,ii) & M.li(k2,i0)~=M.li(k2,ii) );
    if isempty(j)
      continue
    end
    if alg==1 % slow but accurate algorithm, testing ordering constraint on image line segments
      i = [i ii2(ordering_constraint_pairwise(M.li(k,i0),M.li(k,ii2),l(k),P(k)))];
    else % faster but less accurate algorithm, testing ordering constraint on reconstructed 3D line segments
      i = [i ii2(ordering_constraint_pairwise_3d(normx(vgg_wedge(P{k1})),normx(vgg_wedge(P{k2})),M.X(:,i0),M.Y(:,i0),M.X(:,ii2),M.Y(:,ii2)))];
    end
    ii1 = setminus(ii1,i);
  end
  ii = setminus(ii,i);
end

return